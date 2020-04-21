package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.DuplicateSetWalker;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.util.Random;

@CommandLineProgramProperties(
        summary = "Discard a set fraction of duplicate sets from a UMI-grouped ba",
        oneLineSummary = "Discard a set fraction of duplicate sets from a UMI-grouped bam",
        programGroup = ReadDataManipulationProgramGroup.class
)
/**
 * Given a bam grouped by the same unique molecular identifier (UMI), this tool drops a fraction of duplicate sets and returns a new bam.
 * A duplicate set refers to a group of reads whose fragments start at and end at the same coordinate and share the same UMI.
 *
 * The input bam must have been sorted by UMI using FGBio GroupReadsByUmi (http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html).
 *
 * Use this tool to create, for instance, an insilico mixture of duplex-sequenced samples to simulate tumor subclones.
 * Suppose you wish to simulate a tumor sample in which 5% cells share a common set of somatic mutations
 * in addition to ones common to the entire cell population.
 *
 * If you randomly drop 5% of reads in sample A and 95% of reads in sample B and merge the reduced bams,
 * the resulting mixture skews the family-size distribution to the left. Here the family size refers to the
 * number of sequenced duplicate reads that share the same UMI.
 *
 * To see this, take a cancer sample, in which 5% of cells (i.e. a subclone) share a unique set of somatic mutations,
 * that was processed with duplex-UMIs (i.e. UMIs on both adapters) and high rounds of PCR. Suppose we have the sequence-ready
 * libraries of this sample attached to and amplified on the flowcell. Now, sort the flowcell lawn such that the
 * 5% subclone moves near the top of the flowcell. This subclone must have the same family-size distribution as
 * the rest of the flowcell, at about 5% of the library complexity compared to the entire flowcell.
 *
 * Now imagine replacing this subclone with 5% of the *entire* flowcell from another sample prepared and sequenced similarly.
 * The library complexity of these "graft" reads is higher than that of the original, and, consequently, with other parameters
 * such as the number of PCR cycles and sequencing depth fixed, its family distribution would be skewed left---that is, the family size
 * would be smaller than it should be.
 *
 * This tool address the above problem by dropping a set fraction of _duplicate sets_, rather than reads, at random.
 * Implicit in this approach is that a read and its mate are dropped or retained together.
 * While trivial when the input bam is sorted by UMI and query name, this is far from trivial when one attempts
 * to downsample reads naively with a tool like {@link PrintReads}.
 *
 * Example Usage:
 *
 * Keep 95 percent of the molecules.
 *
 * gatk DownsampleByDuplicateSet \ \
 * -I umiGrouped.bam \
 * --fraction-to-keep 0.95 \
 * -O umiGrouped_0.95.bam
 **/
public class DownsampleByDuplicateSet extends DuplicateSetWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "")
    public GATKPathSpecifier outputBam;

    public static final String FRACTION_TO_KEEP_NAME = "fraction-to-keep";
    @Argument(fullName = FRACTION_TO_KEEP_NAME, doc = "This fraction of duplicate sets in the input bam will be retained", minValue = 0.0, maxValue = 1.0)
    public double fractionToKeep;

    private static final int RANDOM_SEED = 142;
    private Random rng;
    private int numDuplicateReadSets;
    private int numReads;
    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        rng = new Random(RANDOM_SEED);
        outputWriter = createSAMWriter(outputBam.toPath(), false);
    }

    @Override
    public void apply(ReadSetWithSharedUMI readSetWithSharedUMI, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (rng.nextDouble() < fractionToKeep){
            readSetWithSharedUMI.getReads().forEach(r -> outputWriter.addRead(r));
            numReads += readSetWithSharedUMI.getReads().size();
            numDuplicateReadSets += 1;
        }
    }

    @Override
    public Object onTraversalSuccess(){
        outputWriter.close();
        logger.info(String.format("Wrote %d reads", numReads));
        logger.info(String.format("Wrote %d duplicate read sets", numDuplicateReadSets));
        return "SUCCESS";
    }

    @Override
    public void closeTool(){
        if (outputWriter != null) {
            outputWriter.close();
        }
    }
}
