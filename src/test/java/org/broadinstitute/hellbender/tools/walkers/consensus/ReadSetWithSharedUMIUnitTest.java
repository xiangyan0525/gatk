package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMTag;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class ReadSetWithSharedUMIUnitTest extends BaseTest {
    // Parameters for the default reads
    final String contig = "3";
    final int numReadPairs = 30;
    final int moleculeNumber = 10;
    final int readLength = 146;

    @Test
    public void testAddRead(){
        final List<GATKRead> reads = makeReads(numReadPairs, contig, moleculeNumber, "read");
        final ReadSetWithSharedUMI set = new ReadSetWithSharedUMI(reads.get(0));
        reads.stream().skip(1).forEach(r -> set.addRead(r));

        Assert.assertEquals(set.getReads().size(), reads.size());
        // Adding a read from a different contig (even if the molecule ID somehow matches) throws an error
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddDifferentMolecule(){
        final List<GATKRead> reads = makeReads(1, contig, moleculeNumber, "one_molecule");
        final ReadSetWithSharedUMI set = new ReadSetWithSharedUMI(reads.get(0));
        reads.stream().skip(1).forEach(set::addRead);

        // Add reads from a different molecule
        final List<GATKRead> readsFromAnotherMolecule = makeReads(1, contig, moleculeNumber + 1, "another_molecule");
        readsFromAnotherMolecule.forEach(set::addRead); // Throws an exception
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddDifferentContig(){
        final List<GATKRead> reads = makeReads(1, contig, moleculeNumber, "one_molecule");
        final ReadSetWithSharedUMI set = new ReadSetWithSharedUMI(reads.get(0));
        reads.stream().skip(1).forEach(set::addRead);

        // Add reads from a different molecule
        final String anotherContig = "4";
        final List<GATKRead> readsFromAnotherMolecule = makeReads(1, anotherContig, moleculeNumber + 1, "another_molecule");
        readsFromAnotherMolecule.forEach(set::addRead); // Throws an exception
    }

    private List<GATKRead> makeReads(final int numReadPairs, final String contig, final int moleculeNumber,
                                     final String readName){
        final List<GATKRead> reads = new ArrayList<>(numReadPairs);
        final List<GATKRead> mates = new ArrayList<>(numReadPairs);

        for (int i = 0; i < numReadPairs; i++){
            // Forward read
            final int readStart = 3000;
            final String readStrand = i % 2 == 0 ? "A" : "B";
            final MoleculeID id = new MoleculeID(moleculeNumber, readStrand);
            final GATKRead read = ArtificialReadUtils.createSamBackedRead(readName + i, contig, readStart, readLength);
            read.setAttribute(SAMTag.MI.name(), id.getSAMField());
            read.setIsReverseStrand(false);

            // Reverse read
            final int insertSize = 1000;
            final int mateStart = readStart + insertSize - readLength;
            final GATKRead mate = ArtificialReadUtils.createSamBackedRead(readName + "_mate" + i, contig, mateStart, readLength);
            final MoleculeID mateId = new MoleculeID(moleculeNumber, readStrand == "A" ? "B" : "A");
            mate.setAttribute(SAMTag.MI.name(), mateId.getSAMField());
            mate.setIsReverseStrand(true);

            if (readStrand.equals("A")){
                read.setIsFirstOfPair(); // F1R2
                mate.setIsSecondOfPair();
            } else {
                read.setIsSecondOfPair(); // F2R1
                mate.setIsFirstOfPair();
            }

            reads.add(read);
            mates.add(mate);
        }

        final List<GATKRead> results = new ArrayList<>(numReadPairs*2);
        results.addAll(reads);
        results.addAll(mates); // Add reads and mates in order so that the resulting list is ordered by coordinate
        return results;
    }
}