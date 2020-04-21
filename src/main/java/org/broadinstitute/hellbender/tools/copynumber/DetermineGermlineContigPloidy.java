package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineContigPloidyHybridADVIArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineContigPloidyModelArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CoveragePerContigCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CoveragePerContig;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils.streamOfSubsettedAndValidatedReadCounts;

/**
 * Determines the integer ploidy state of all contigs for germline samples given counts data. These should be either
 * HDF5 or TSV count files generated by {@link CollectReadCounts}; TSV files may be compressed (e.g., with bgzip),
 * but must then have filenames ending with the extension .gz.  See the documentation for the {@code input} argument
 * for details on enabling streaming of indexed count files from Google Cloud Storage.
 *
 * <h3>Introduction</h3>
 *
 * <p>Germline karyotyping is a frequently performed task in bioinformatics pipelines, e.g. for sex determination and
 * aneuploidy identification. This tool uses counts data for germline karyotyping. </p>
 *
 * <p>Performing germline karyotyping using counts data requires calibrating ("modeling") the technical coverage bias
 * and variance for each contig. The Bayesian model and the associated inference scheme implemented in
 * {@link DetermineGermlineContigPloidy} includes provisions for inferring and explaining away much of the technical
 * variation. Furthermore, karyotyping confidence is automatically adjusted for individual samples and contigs.</p>
 *
 * <p>Running {@link DetermineGermlineContigPloidy} is the first computational step in the GATK germline CNV calling
 * pipeline. It provides a baseline ("default") copy-number state for each contig/sample with respect to which the
 * probability of alternative states is allocated.</p>
 *
 * <h3>Python environment setup</h3>
 *
 * <p>The computation done by this tool, aside from input data parsing and validation, is performed outside of the Java
 * Virtual Machine and using the <em>gCNV computational python module</em>, namely {@code gcnvkernel}. It is crucial that
 * the user has properly set up a python conda environment with {@code gcnvkernel} and its dependencies
 * installed. If the user intends to run {@link DetermineGermlineContigPloidy} using one of the official GATK Docker images,
 * the python environment is already set up. Otherwise, the environment must be created and activated as described in the
 * main GATK README.md file.</p>
 *
 * <p>Advanced users may wish to set the <code>THEANO_FLAGS</code> environment variable to override the GATK theano
 * configuration. For example, by running
 * <code>THEANO_FLAGS="base_compiledir=PATH/TO/BASE_COMPILEDIR" gatk DetermineGermlineContigPloidy ...</code>, users can specify
 * the theano compilation directory (which is set to <code>$HOME/.theano</code> by default).  See theano documentation
 * at <a href="http://deeplearning.net/software/theano/library/config.html">
 *     http://deeplearning.net/software/theano/library/config.html</a>.
 * </p>
 *
 * <h3>Tool run modes</h3>
 *
 * <p>This tool has two operation modes as described below:</p>
 * <dl>
 *     <dt>COHORT mode:</dt>
 *
 *     <dd>If a ploidy model parameter path is not provided via the {@code model} argument, the tool will run in
 *     the COHORT mode. In this mode, ploidy model parameters (e.g. coverage bias and variance for each contig) are
 *     inferred, along with baseline contig ploidy states of each sample. It is possible to run the tool over a subset
 *     of all intervals present in the input count files, which can be specified by -L; this can be used to pass a
 *     filtered interval list produced by {@link FilterIntervals} to mask intervals from modeling. Intervals may also be
 *     blacklisted using -XL. The specified intervals that result from resolving -L/-XL inputs must be exactly present
 *     in all of the input count files.
 *
 *     <p>A TSV file specifying prior probabilities for each integer ploidy state and for each contig is required in this
 *     mode and must be specified via the {@code contig-ploidy-priors} argument. The following shows an example of
 *     such a table:</p>
 *     <br>
 *     <table border="1" width="80%">
 *         <tr>
 *             <td>CONTIG_NAME</td> <td>PLOIDY_PRIOR_0</td> <td>PLOIDY_PRIOR_1</td> <td>PLOIDY_PRIOR_2</td>
 *             <td>PLOIDY_PRIOR_3</td>
 *         </tr>
 *         <tr>
 *            <td>1</td> <td>0.01</td> <td>0.01</td> <td>0.97</td> <td>0.01</td>
 *         </tr>
 *         <tr>
 *            <td>2</td> <td>0.01</td> <td>0.01</td> <td>0.97</td> <td>0.01</td>
 *         </tr>
 *         <tr>
 *             <td>X</td> <td>0.01</td> <td>0.49</td> <td>0.49</td> <td>0.01</td>
 *         </tr>
 *         <tr>
 *            <td>Y</td> <td>0.50</td> <td>0.50</td> <td>0.00</td> <td>0.00</td>
 *         </tr>
 *     </table>
 *     <p>Note that the contig names appearing under {@code CONTIG_NAME} column must match contig names in the input
 *     counts files, and all contigs appearing in the input counts files must have a corresponding entry in the priors
 *     table. The order of contigs is immaterial in the priors table. The highest ploidy state is determined by the
 *     prior table (3 in the above example). A ploidy state can be strictly forbidden by setting its prior probability
 *     to 0. For example, the Y contig in the above example can only assume 0 and 1 ploidy states.</p>
 *
 *     <p>The tool output in the COHORT mode will contain two subdirectories, one ending with "-model" and the other
 *     ending with "-calls". The model subdirectory contains the inferred parameters of the ploidy model, which may
 *     be used later on for karyotyping one or more similarly-sequenced samples (see below).
 *
 *     The calls subdirectory contains one subdirectory for each sample, listing various sample-specific
 *     quantities such as the global read-depth, average ploidy, per-contig baseline ploidies, and per-contig
 *     coverage variance estimates.</p></dd>
 *
 *     <dt>CASE mode:</dt>
 *     <dd>If a path containing previously inferred ploidy model parameters is provided via the
 *     {@code model} argument, then the tool will run in the CASE mode. In this mode, the parameters of the ploidy
 *     model are loaded from the provided directory and only sample-specific quantities are inferred. The modeled
 *     intervals are then specified by a file contained in the model directory, all interval-related arguments are
 *     ignored in this mode, and all model intervals must be present in all of the input count files. The tool output
 *     in the CASE mode is only the "-calls" subdirectory and is organized similarly to the COHORT mode.
 *
 *      <p>In the CASE mode, the contig ploidy prior table is taken directly from the provided model parameters
 *      path and must be not provided again.</p></dd>
 * </dl>
 *
 * <h3>Important Remarks</h3>
 * <dl>
 *     <dt>Choice of hyperparameters:</dt>
 *     <dd><p>The quality of ploidy model parametrization and the sensitivity/precision of germline karyotyping are
 *     sensitive to the choice of model hyperparameters, including standard deviation of mean contig coverage bias
 *     (set using the {@code mean-bias-standard-deviation} argument), mapping error rate
 *     (set using the {@code mapping-error-rate} argument), and the typical scale of contig- and sample-specific
 *     unexplained variance (set using the {@code global-psi-scale} and {@code sample-psi-scale} arguments,
 *     respectively). It is crucial to note that these hyperparameters are <em>not</em> universal
 *     and must be tuned for each sequencing protocol and properly set at runtime.</p></dd>
 *
 *     <dt>Mosaicism and fractional ploidies:</dt>
 *     <dd><p>The model underlying this tool assumes integer ploidy states (in contrast to fractional/variable ploidy
 *     states). Therefore, it is to be used strictly on germline samples and for the purpose of sex determination,
 *     autosomal aneuploidy detection, or as a part of the GATK germline CNV calling pipeline. The presence of large somatic
 *     events and mosaicism (e.g., sex chromosome loss and somatic trisomy) will naturally lead to unreliable
 *     results. We strongly recommended inspecting genotyping qualities (GQ) from the tool output and considering to drop
 *     low-GQ contigs in downstream analyses. Finally, given the Bayesian status of this tool, we suggest including as many
 *     high-quality germline samples as possible for ploidy model parametrizaton in the COHORT mode. This will downplay
 *     the role of questionable samples and will yield a more reliable estimation of genuine sequencing biases.</p></dd>
 *
 *     <dt>Coverage-based germline karyotyping:</dt>
 *     <dd><p>Accurate germline karyotyping requires incorporating SNP allele-fraction data and counts data in a
 *     unified probabilistic model and is beyond the scope of the present tool. The current implementation only uses
 *     counts data for karyotyping and while being fast, it may not provide the most reliable results.</p></dd>
 *
 * <h3>Usage examples</h3>
 *
 * <p>COHORT mode:</p>
 * <pre>
 * gatk DetermineGermlineContigPloidy \
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --contig-ploidy-priors a_valid_ploidy_priors_table.tsv
 *   --output output_dir \
 *   --output-prefix normal_cohort
 * </pre>
 *
 * <p>COHORT mode (with optional interval filtering):</p>
 * <pre>
 * gatk DetermineGermlineContigPloidy \
 *   -L intervals.interval_list \
 *   -XL blacklist_intervals.interval_list \
 *   --interval-merging-rule OVERLAPPING_ONLY \
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --contig-ploidy-priors a_valid_ploidy_priors_table.tsv
 *   --output output_dir \
 *   --output-prefix normal_cohort
 * </pre>
 *
 * <p>CASE mode:</p>
 * <pre>
 * gatk DetermineGermlineContigPloidy \
 *   --model a_valid_ploidy_model_dir
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --output output_dir \
 *   --output-prefix normal_case
 * </pre>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Determines the baseline contig ploidy for germline samples given counts data",
        oneLineSummary = "Determines the baseline contig ploidy for germline samples given counts data",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class DetermineGermlineContigPloidy extends CommandLineProgram {
    public enum RunMode {
        COHORT, CASE
    }

    public static final String COHORT_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT = "cohort_determine_ploidy_and_depth.py";
    public static final String CASE_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT = "case_determine_ploidy_and_depth.py";

    //name of the interval file output by the python code in the model directory
    public static final String INPUT_MODEL_INTERVAL_FILE = "interval_list.tsv";

    public static final String MODEL_PATH_SUFFIX = "-model";
    public static final String CALLS_PATH_SUFFIX = "-calls";

    public static final String CONTIG_PLOIDY_PRIORS_FILE_LONG_NAME = "contig-ploidy-priors";

    @Argument(
            doc = "Input paths for read-count files containing integer read counts in genomic intervals for all samples.  " +
                    "All intervals specified via -L/-XL must be contained; " +
                    "if none are specified, then intervals must be identical and in the same order for all samples.  " +
                    "If read-count files are given by Google Cloud Storage paths, " +
                    "have the extension .counts.tsv or .counts.tsv.gz, " +
                    "and have been indexed by IndexFeatureFile, " +
                    "only the specified intervals will be queried and streamed; " +
                    "this can reduce disk usage by avoiding the complete localization of all read-count files.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            minElements = 1
    )
    private List<String> inputReadCountPaths = new ArrayList<>();

    @Argument(
            doc = "Input file specifying contig-ploidy priors.  If only a single sample is specified, this input should not be provided.  " +
                    "If multiple samples are specified, this input is required.",
            fullName = CONTIG_PLOIDY_PRIORS_FILE_LONG_NAME,
            optional = true
    )
    private File inputContigPloidyPriorsFile;

    @Argument(
            doc = "Input ploidy-model directory.  If only a single sample is specified, this input is required.  " +
                    "If multiple samples are specified, this input should not be provided.",
            fullName = CopyNumberStandardArgument.MODEL_LONG_NAME,
            optional = true
    )
    private File inputModelDir;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.  This will be created if it does not exist.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputDir;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection
            = new OptionalIntervalArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineContigPloidyModelArgumentCollection germlineContigPloidyModelArgumentCollection =
            new GermlineContigPloidyModelArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineContigPloidyHybridADVIArgumentCollection germlineContigPloidyHybridADVIArgumentCollection =
            new GermlineContigPloidyHybridADVIArgumentCollection();

    private RunMode runMode;
    private SimpleIntervalCollection specifiedIntervals;
    private File specifiedIntervalsFile;

    @Override
    protected void onStartup() {
        /* check for successful import of gcnvkernel */
        PythonScriptExecutor.checkPythonEnvironmentForPackage("gcnvkernel");
    }

    @Override
    protected Object doWork() {
        validateArguments();

        setModeAndResolveIntervals();

        //read in count files and output intervals and samples x coverage-per-contig table to temporary files
        specifiedIntervalsFile = IOUtils.createTempFile("intervals", ".tsv");
        specifiedIntervals.write(specifiedIntervalsFile);
        final File samplesByCoveragePerContigFile = IOUtils.createTempFile("samples-by-coverage-per-contig", ".tsv");
        writeSamplesByCoveragePerContig(samplesByCoveragePerContigFile);

        //call python inference code
        final boolean pythonReturnCode = executeDeterminePloidyAndDepthPythonScript(
                samplesByCoveragePerContigFile, specifiedIntervalsFile);

        if (!pythonReturnCode) {
            throw new UserException("Python return code was non-zero.");
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void validateArguments() {
        germlineContigPloidyModelArgumentCollection.validate();
        germlineContigPloidyHybridADVIArgumentCollection.validate();

        Utils.validateArg(inputReadCountPaths.size() == new HashSet<>(inputReadCountPaths).size(),
                "List of input read-count files cannot contain duplicates.");

        inputReadCountPaths.forEach(CopyNumberArgumentValidationUtils::validateInputs);
        CopyNumberArgumentValidationUtils.validateInputs(
                inputContigPloidyPriorsFile,
                inputModelDir);
        CopyNumberArgumentValidationUtils.checkForSingletonIntervals(specifiedIntervals);
        Utils.nonEmpty(outputPrefix);
        CopyNumberArgumentValidationUtils.validateAndPrepareOutputDirectories(outputDir);
    }

    private void setModeAndResolveIntervals() {
        if (inputModelDir != null) {
            runMode = RunMode.CASE;
            logger.info("A contig-ploidy model was provided, running in case mode...");
            if (inputContigPloidyPriorsFile != null) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in case mode, " +
                        "but contig-ploidy priors were provided.");
            }
            if (intervalArgumentCollection.intervalsSpecified()) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in CASE mode, " +
                        "but intervals were provided.");
            }
            //intervals are retrieved from the input model directory
            specifiedIntervalsFile = new File(inputModelDir, INPUT_MODEL_INTERVAL_FILE);
            CopyNumberArgumentValidationUtils.validateInputs(specifiedIntervalsFile);
            specifiedIntervals = new SimpleIntervalCollection(specifiedIntervalsFile);
        } else {
            runMode = RunMode.COHORT;
            logger.info("No contig-ploidy model was provided, running in cohort mode...");
            Utils.validateArg(inputReadCountPaths.size() > 1, "At least two samples must be provided in " +
                    "COHORT mode.");
            if (inputContigPloidyPriorsFile == null){
                throw new UserException.BadInput("Contig-ploidy priors must be provided in cohort mode.");
            }

            //get sequence dictionary and intervals from the first read-count file to use to validate remaining files
            //(this first file is read again below, which is slightly inefficient but is probably not worth the extra code)
            final String firstReadCountPath = inputReadCountPaths.get(0);
            specifiedIntervals = CopyNumberArgumentValidationUtils.resolveIntervals(
                    firstReadCountPath, intervalArgumentCollection, logger);
        }
    }

    private void writeSamplesByCoveragePerContig(final File samplesByCoveragePerContigFile) {
        //calculate coverage per contig and construct record for each sample
        logger.info("Validating and aggregating coverage per contig from input read-count files...");
        final List<String> contigs = specifiedIntervals.getRecords().stream()
                .map(SimpleInterval::getContig)
                .distinct()
                .collect(Collectors.toList());
        final List<CoveragePerContig> coveragePerContigs =
                streamOfSubsettedAndValidatedReadCounts(inputReadCountPaths, specifiedIntervals, logger)
                        .map(subsetReadCounts ->
                                new CoveragePerContig(
                                        subsetReadCounts.getMetadata().getSampleName(),
                                        subsetReadCounts.getRecords().stream()
                                                .collect(Collectors.groupingBy(
                                                        SimpleCount::getContig,
                                                        LinkedHashMap::new,
                                                        Collectors.summingInt(SimpleCount::getCount)))))
                        .collect(Collectors.toList());
        new CoveragePerContigCollection(specifiedIntervals.getMetadata(), coveragePerContigs, contigs)
                .write(samplesByCoveragePerContigFile);
    }

    private boolean executeDeterminePloidyAndDepthPythonScript(final File samplesByCoveragePerContigFile,
                                                               final File intervalsFile) {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final String outputDirArg = CopyNumberArgumentValidationUtils.addTrailingSlashIfNecessary(outputDir.getAbsolutePath());
        //note that the samples x coverage-by-contig table is referred to as "metadata" by gcnvkernel
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                "--sample_coverage_metadata=" + CopyNumberArgumentValidationUtils.getCanonicalPath(samplesByCoveragePerContigFile),
                "--output_calls_path=" + CopyNumberArgumentValidationUtils.getCanonicalPath(outputDirArg + outputPrefix + CALLS_PATH_SUFFIX)));
        arguments.addAll(germlineContigPloidyModelArgumentCollection.generatePythonArguments(runMode));
        arguments.addAll(germlineContigPloidyHybridADVIArgumentCollection.generatePythonArguments());

        final String script;
        if (runMode == RunMode.COHORT) {
            script = COHORT_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT;
            arguments.add("--interval_list=" + CopyNumberArgumentValidationUtils.getCanonicalPath(intervalsFile));
            arguments.add("--contig_ploidy_prior_table=" + CopyNumberArgumentValidationUtils.getCanonicalPath(inputContigPloidyPriorsFile));
            arguments.add("--output_model_path=" + CopyNumberArgumentValidationUtils.getCanonicalPath(outputDirArg + outputPrefix + MODEL_PATH_SUFFIX));
        } else {
            script = CASE_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT;
            arguments.add("--input_model_path=" + CopyNumberArgumentValidationUtils.getCanonicalPath(inputModelDir));
        }
        return executor.executeScript(
                new Resource(script, GermlineCNVCaller.class),
                null,
                arguments);
    }
}
