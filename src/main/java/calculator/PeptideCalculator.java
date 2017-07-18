package calculator;

import model.SecondaryStructureType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;

import java.util.Set;

public class PeptideCalculator {

    public void structureToPeptide(Structure structure, int peptideBaseSize, int peptideMinSize,
                                   int peptideMaxSize, int peptideOverlapBase, int peptideOverlapMax,
                                   Set<SecondaryStructureType> secondaryStructureTypes, boolean secondaryStructureSensitivity) throws StructureException {
        SecStrucCalc secondaryStructureCalculator = new SecStrucCalc();
        secondaryStructureCalculator.calculate(structure ,true);

//        for (structure.chain)
    }
}
