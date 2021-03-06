import calculator.PeptideCalculator;
import model.Peptide;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.secstruc.*;


import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Main {

    public static void main(String[] args) throws IOException,
            StructureException {

        String pdbID = "5A0Y";

        // Only change needed to the DEFAULT Structure loading
        FileParsingParameters params = new FileParsingParameters();
        params.setParseSecStruc(true);

        AtomCache cache = new AtomCache();
        cache.setFileParsingParams(params);

        // Use PDB format, because SS cannot be parsed from mmCIF yet
        cache.setUseMmCif(false);

        // The loaded Structure contains the SS assigned by Author (simple)
        Structure s = cache.getStructure(pdbID);
        PeptideCalculator calculator = new PeptideCalculator();
        Set<SecStrucType> secStrucTypes = new HashSet<SecStrucType>();
        secStrucTypes.add(SecStrucType.helix3);
        secStrucTypes.add(SecStrucType.helix4);
        secStrucTypes.add(SecStrucType.helix5);
        secStrucTypes.add(SecStrucType.bridge);
        List<Map<String, List<Peptide>>> models = calculator.structureToPeptides(s, 10, 8, 12, 2, 4, secStrucTypes, true);

        // Print the Author's assignment (from PDB file)
//        System.out.println("Author's assignment: ");
//        printSecStruc(s);

//        // If the more detailed DSSP prediction is required call this
//        DSSPParser.fetch(pdbID, s, true);
//
//        // Print the assignment residue by residue
//        System.out.println("DSSP assignment: ");
//        printSecStruc(s);

        // finally use BioJava's built in DSSP-like secondary structure assigner
        SecStrucCalc secStrucCalc = new SecStrucCalc();

        // calculate and assign
        List<SecStrucState> l = secStrucCalc.calculate(s,true);
        printSecStruc(s);

    }

    public static void printSecStruc(Structure s){
        List<SecStrucInfo> ssi = SecStrucTools.getSecStrucInfo(s);
        for (SecStrucInfo ss : ssi) {
            System.out.println(ss.getGroup().getChain().getChainID() + " "
                    + ss.getGroup().getResidueNumber() + " "
                    + ss.getGroup().getPDBName() + " -> " + ss.toString());
        }
    }
}
