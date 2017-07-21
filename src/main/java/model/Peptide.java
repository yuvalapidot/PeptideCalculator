package model;

import org.biojava.nbio.structure.Group;

import java.util.List;

/**
 * Created by yuvalapidot.
 */
public class Peptide {

    private List<Group> aminoAcids;
    private int startIndex;
    private int endIndex;
    private int size;

    public Peptide(List<Group> aminoAcids, int startIndex, int size) {
        this.aminoAcids = aminoAcids;
        this.startIndex = startIndex;
        this.size = size;
        this.endIndex = startIndex + size - 1;
    }

    public List<Group> getAminoAcids() {
        return aminoAcids;
    }

    public int getStartIndex() {
        return startIndex;
    }

    public int getEndIndex() {
        return endIndex;
    }

    public int getSize() {
        return size;
    }
}
