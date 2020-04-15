package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.walkers.mutect.UMI;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A container class for a set of reads that share the same unique molecular identifier, as determined by
 * FGBio GroupReadsByUmi (http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)
 */
public class ReadSetWithSharedUMI implements Locatable {
    public static final String FGBIO_MI_TAG_DELIMITER = "/";
    private static final int EMPTY_MOLECULE_ID = -1;
    private int moleculeId = EMPTY_MOLECULE_ID; // TODO: extract an ID class.

    private String contig;
    private int fragmentStart = -1;
    private int fragmentEnd = -1;

    private List<GATKRead> reads;

    public ReadSetWithSharedUMI(){
        reads = new ArrayList<>();
    }

    public ReadSetWithSharedUMI(final GATKRead read){
        this();
        init(read);
        reads.add(read);
    }

    private void init(GATKRead read){
        Utils.validate(moleculeId == -1 || moleculeId == getMoleculeID(read),
                String.format("Inconsistent molecule IDs: Duplicate set id = %s, read molecule id = %s", moleculeId, getMoleculeID(read)));
        setMoleduleId(read);
        contig = read.getContig();
        fragmentStart = read.getStart();
        fragmentEnd = read.getEnd(); // TODO: does this include softclips?
    }

    public List<GATKRead> getReads(){
        return Collections.unmodifiableList(reads);
    }

    public boolean sameMolecule(final GATKRead read){
        return getMoleculeID(read) == moleculeId;
    }

    public void setMoleduleId(GATKRead read){
        moleculeId = getMoleculeID(read);
    }

    public void setMoleduleId(int id){
        moleculeId = id;
    }

    /**
     * Some examples of molecule IDs (MI tag):
     *
     * "0/A" (The first molecule in the bam, top (A) strand)
     * "0/B" (The first molecule in the bam, bottom (B) strand)
     * "99/A" (100th molecule in the bam, top top (A) strand)
     *
     * Top strand is synonymous to "F1R2"
     * Bottom strand is synonymous to "F2R1"
     *
     * Thus only the integer component is relevant for identifying reads that originated from the same molecule.
     * Should the need arise, we could extend this to distinguish between different strands of the same molecule.
     */
    public static int getMoleculeID(final GATKRead read) {
        Utils.validateArg(read.hasAttribute(SAMTag.MI.name()),
                "Reads must have the molecular ID Tag: " + SAMTag.MI.name() + ". Read = " + read.toString());
        final String MITag = read.getAttributeAsString(SAMTag.MI.name());
        return Integer.parseInt(MITag.split(FGBIO_MI_TAG_DELIMITER)[0]);
    }

    /** Returns true if the read was properly added to the duplicate set **/
    public boolean addRead(final GATKRead read){
        if (reads.isEmpty()){
            init(read);
            reads.add(read);
            return true;
        }

        if (sameMolecule(read)){
            reads.add(read);
            if (read.getStart() < fragmentStart){
                fragmentStart = read.getStart();
            }
            if (read.getEnd() > fragmentEnd){
                fragmentEnd = read.getEnd();
            }
            return true;
        } else {
            return false;
        }
    }

    public boolean isEmpty(){
        return reads.isEmpty();
    }

    public int getMoleculeId() {
        if (this.isEmpty()){
            return EMPTY_MOLECULE_ID;
        }
        return moleculeId;
    }

    public SimpleInterval getDuplicateSetInterval(){
        return new SimpleInterval(contig, fragmentStart, fragmentEnd);
    }

    public boolean hasValidInterval(){
        return SimpleInterval.isValid(contig, fragmentStart, fragmentEnd);
    }

    public static List<String> getMolecularIDs(final List<GATKRead> reads) {
        return reads.stream().map(r -> r.getAttributeAsString(ReadSetWithSharedUMI.SAMTag.MI.name()))
                .distinct().collect(Collectors.toList());
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return fragmentStart;
    }

    @Override
    public int getEnd() {
        return fragmentEnd;
    }

    /**
     * A container class for the molecular ID, which consists of an integer ID and a binary strand.
     * For example, Reads with the tags 12/A and 12/B originated from the same DNA fragment before PCR,
     * (i.e. from the same library) but they originated from different strands in that library.
     * Thus one read is F1R2 and the other F2R1.
     */
    private static class MoleculeIDTag {
        private int molecularID;
        private String strand;

        private MoleculeIDTag(final GATKRead read){
            final String MITag = read.getAttributeAsString(SAMTag.MI.name());

            this.molecularID = Integer.parseInt(MITag.split(FGBIO_MI_TAG_DELIMITER)[0]);
            this.strand = MITag.split(FGBIO_MI_TAG_DELIMITER)[1];
        }

        public int getMoleculeID() {
            return molecularID;
        }

        public String getStrand() {
            return strand;
        }
    }
}
