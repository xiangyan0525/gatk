package org.broadinstitute.hellbender.engine;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.consensus.MoleculeID;
import org.broadinstitute.hellbender.tools.walkers.consensus.ReadSetWithSharedUMI;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A walker that processes duplicate reads that share the same UMI as a single unit.
 *
 * This tool assumes that the input bam has been sorted by UMI with FGBio GroupReadsByUmi:
 * http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html
 */
public abstract class DuplicateSetWalker extends ReadWalker {
    public static final String MIN_REQUIRED_READS_NAME = "min-reads";
    public static final String MIN_REQUIRED_READS_PER_STRAND_NAME = "min-per-strand-reads";

    private static final int DEFAULT_MINIMUM_READS_PER_SET = 1;
    private static final int DEFAULT_MINIMUM_READS_PER_STRAND = 0;

    @Argument(fullName = MIN_REQUIRED_READS_NAME, doc = "The mininum total number of reads required in the set", optional = true)
    private int minimumRequiredReadsPerUMI = DEFAULT_MINIMUM_READS_PER_SET;

    /** The user may keep only read sets with both strands by setting this argument to a positive number **/
    @Argument(fullName = MIN_REQUIRED_READS_PER_STRAND_NAME, doc = "The mininum total number of reads in each strand", optional = true)
    private int minimumRequiredReadsPerStrand = DEFAULT_MINIMUM_READS_PER_STRAND;

    private ReadSetWithSharedUMI currentReadSetWithSharedUMI = null;

    /***
     * FGBio GroupByUMI returns reads sorted by molecular ID: For example, the input bam may look like
     * read1: ... MI:Z:0/A ...
     * read2: ... MI:Z:0/A ...
     * read3: ... MI:Z:0/B ...
     * read4: ... MI:Z:0/B ...
     * read5: ... MI:Z:1/A ...
     * read6: ... MI:Z:1/B ...
     * read7: ... MI:Z:1/B ...
     *
     * Thus it's sufficient to go through the reads in order and collect them in a list until
     * we encounter the next molecular ID, at which point we pass the list to the {@code apply}
     * method of the child class and clear the {@code currentDuplicateSet} variable.
     *
     * Marked final to discourage subclassing. A subclass must override the other apply() method that takes in the DuplicateSet.
     * Marked public to match the parent method signature.
     */
    @Override
    public final void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (currentReadSetWithSharedUMI == null){ // evaluate to true for the very first read
            currentReadSetWithSharedUMI = new ReadSetWithSharedUMI(read);
            return;
        }

        final int readMoleculeNumber = MoleculeID.getMoleculeNumberOfRead(read);
        final int duplicateSetMoleculeNumber = currentReadSetWithSharedUMI.getMoleculeNumber();

        // If the incoming read has the molecular id less than that of the currentDuplicateSet,
        // the bam is not sorted properly by the MI tag
        if (duplicateSetMoleculeNumber > readMoleculeNumber){
            throw new UserException("The input bam must be sorted by the molecularID tag.");
        }

        // Check for greater than 0 because before apply() is called on the first read, the set would be empty and we use -1 as a placeholder
        if (duplicateSetMoleculeNumber < readMoleculeNumber) {
            // If the currentDuplicateSet is empty or does not match the incoming read, we've reached the end of the current set.
            // Call the apply() method on the current set and start a new set.
            if (rejectSet(currentReadSetWithSharedUMI)){
                currentReadSetWithSharedUMI = new ReadSetWithSharedUMI(read);
                return;
            }

            apply(currentReadSetWithSharedUMI,
                    new ReferenceContext(reference, currentReadSetWithSharedUMI.getInterval()), // Will create an empty ReferenceContext if reference or readInterval == null
                    new FeatureContext(features, currentReadSetWithSharedUMI.getInterval()));
            currentReadSetWithSharedUMI = new ReadSetWithSharedUMI(read);
            return;
        }


        // the incoming read has the same UMI as the current set
        currentReadSetWithSharedUMI.addRead(read);
    }

    /**
     * A subclass must specify how to process the duplicate sets by overriding this method.
     *
     * @param readSetWithSharedUMI A set of reads with the matching UMIs with the same fragment start and end
     * @param referenceContext A reference context object over the intervals determined by the duplicate set.
     * @param featureContext Entries from a secondary feature file (e.g. vcf) if provided
     *
     */
    public abstract void apply(ReadSetWithSharedUMI readSetWithSharedUMI, ReferenceContext referenceContext, FeatureContext featureContext );


    /**
     * Returns true for duplicate sets that does not meet required criteria for further processing.
     *
     * We encourage the user override this method to meet their needs.
     */
    protected boolean rejectSet(final ReadSetWithSharedUMI readSetWithSharedUMI){
        // Check that the set contains the minimum required number of reads in each strand
        final Pair<Integer, Integer> strandCounts = MoleculeID.countStrands(readSetWithSharedUMI.getReads());
        if (Math.min(strandCounts.getLeft(), strandCounts.getRight()) < minimumRequiredReadsPerStrand){
            return true;
        }

        // Check that the read set is paired (checking for a sufficient condition)
        if (readSetWithSharedUMI.getReads().size() % 2 == 1){
            return true;
        }

        // Check that the total number of reads (from both strands) exceeds a specified threshold
        return readSetWithSharedUMI.getReads().size() < minimumRequiredReadsPerUMI;
    }

    @Override
    public void postProcess(){
        if (currentReadSetWithSharedUMI.getReads().size() > 0){
            apply(currentReadSetWithSharedUMI,
                    new ReferenceContext(reference, currentReadSetWithSharedUMI.getInterval()),
                    new FeatureContext(features, currentReadSetWithSharedUMI.getInterval()));
        }
    }
}
