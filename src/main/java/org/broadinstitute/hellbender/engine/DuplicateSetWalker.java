package org.broadinstitute.hellbender.engine;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.consensus.ReadSetWithSharedUMI;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A walker that processes duplicate reads that share the same UMI as a single unit.
 *
 * This tool assumes that the input bam has been sorted by UMI with FGBio GroupReadsByUmi:
 * http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html
 */
public abstract class DuplicateSetWalker extends ReadWalker {
    public static final String MIN_REQUIRED_READS_NAME = "min";
    private static final int DEFAULT_MINIMUM_READS_IN_DUPLICATE_SET = 1;

    @Argument(fullName = MIN_REQUIRED_READS_NAME, optional = true)
    private int minimumRequiredReadsInDuplicateSet = DEFAULT_MINIMUM_READS_IN_DUPLICATE_SET;

    private ReadSetWithSharedUMI currentReadSetWithSharedUMI = new ReadSetWithSharedUMI();

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
        final int readMoleculeID = ReadSetWithSharedUMI.getMoleculeID(read);
        final int duplicateSetMoleculeId = currentReadSetWithSharedUMI.getMoleculeId();

        // If the incoming read has the molecular id less than that of the currentDuplicateSet,
        // the bam is not sorted properly by the MI tag
        if (duplicateSetMoleculeId > readMoleculeID){
            throw new UserException("The input bam must be sorted by the molecularID tag.");
        } else if (duplicateSetMoleculeId < readMoleculeID) {
            // If the currentDuplicateSet is empty or does not match the incoming read,
            // call apply() method on the current duplicate set, and replace it with a new set with the most recent read in it
            if (rejectDuplicateSet(currentReadSetWithSharedUMI)){
                currentReadSetWithSharedUMI = new ReadSetWithSharedUMI(read);
                return;
            }

            apply(currentReadSetWithSharedUMI,
                    new ReferenceContext(reference, currentReadSetWithSharedUMI.getDuplicateSetInterval()), // Will create an empty ReferenceContext if reference or readInterval == null
                    new FeatureContext(features, currentReadSetWithSharedUMI.getDuplicateSetInterval()));
            currentReadSetWithSharedUMI = new ReadSetWithSharedUMI(read);
        } else {
            // the molecule ID of the read matches that of the current duplicate set
            Utils.validate(currentReadSetWithSharedUMI.addRead(read), "Adding a read that doesn't have a matching molecular ID tag");
        }
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
     * A subclass may override this method to meet its needs.
     */
    protected boolean rejectDuplicateSet(final ReadSetWithSharedUMI readSetWithSharedUMI){
        return !readSetWithSharedUMI.hasValidInterval() || readSetWithSharedUMI.getReads().size() < minimumRequiredReadsInDuplicateSet;
    }

    @Override
    public void postProcess(){
        if (currentReadSetWithSharedUMI.getReads().size() > 0){
            apply(currentReadSetWithSharedUMI,
                    new ReferenceContext(reference, currentReadSetWithSharedUMI.getDuplicateSetInterval()),
                    new FeatureContext(features, currentReadSetWithSharedUMI.getDuplicateSetInterval()));
        }
    }
}
