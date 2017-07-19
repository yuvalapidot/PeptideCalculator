package calculator;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.secstruc.SecStrucState;
import org.biojava.nbio.structure.secstruc.SecStrucType;

import java.util.*;

public class PeptideCalculator {

    private static final String SECSTRUCT_PROPERTY_KEY = "secstruc";

    public List<Map<String, List<List<Group>>>> structureToPeptides(Structure structure, int peptideBaseSize,
                                                              int peptideMinSize, int peptideMaxSize,
                                                              int peptideOverlapBase, int peptideOverlapMax,
                                                              Set<SecStrucType> secondaryStructureTypes,
                                                              boolean secondaryStructureSensitivity) throws StructureException {
        assignSecStructureProperty(structure);
        List<Map<String, List<List<Group>>>> models = new ArrayList<Map<String, List<List<Group>>>>();
        for (int modelIndex = 0; modelIndex < structure.nrModels(); modelIndex++) {
            List<Chain> model = structure.getModel(modelIndex);
            models.add(modelToPeptides(model, peptideBaseSize, peptideMinSize,
            peptideMaxSize, peptideOverlapBase, peptideOverlapMax,
            secondaryStructureTypes, secondaryStructureSensitivity));
        }
        return models;
    }

    private void assignSecStructureProperty(Structure structure) throws StructureException {
        SecStrucCalc secStrucCalc = new SecStrucCalc();
        secStrucCalc.calculate(structure,true);
    }

    private Map<String, List<List<Group>>> modelToPeptides(List<Chain> model, int peptideBaseSize, int peptideMinSize, int peptideMaxSize, int peptideOverlapBase, int peptideOverlapMax, Set<SecStrucType> secondaryStructureTypes, boolean secondaryStructureSensitivity) {
        Map<String, List<List<Group>>> chainMap = new HashMap<String, List<List<Group>>>();
        for (Chain chain : model) {
            List<List<Group>> peptides = chainToPeptides(chain, peptideBaseSize, peptideMinSize,
                    peptideMaxSize, peptideOverlapBase, peptideOverlapMax,
                    secondaryStructureTypes, secondaryStructureSensitivity);
            chainMap.put(chain.getChainID(), peptides);
        }
        return chainMap;
    }

    private List<List<Group>> chainToPeptides(Chain chain, int peptideBaseSize, int peptideMinSize, int peptideMaxSize, int peptideOverlapBase, int peptideOverlapMax, Set<SecStrucType> secondaryStructureTypes, boolean secondaryStructureSensitivity) {
        List<List<Group>> peptides = new ArrayList<List<Group>>();
        char[] characterArray = formatChainAsCharacterArray(chain, secondaryStructureTypes, secondaryStructureSensitivity);
        return peptides;
    }

    private char[] formatChainAsCharacterArray(Chain chain, Set<SecStrucType> secondaryStructureTypes, boolean secondaryStructureSensitivity) {
        char[] characterArray = new char[chain.getSeqResLength()];
        int groupIndex = 0;
        for (Group group : chain.getSeqResGroups()) {
            Object propertyValue = group.getProperty(SECSTRUCT_PROPERTY_KEY);
            if (propertyValue == null || !(propertyValue instanceof SecStrucState)) {
                characterArray[groupIndex] = 0;
            } else {
                SecStrucState secStrucState = (SecStrucState) propertyValue;
                if (secondaryStructureTypes.contains(secStrucState.getType())) {
                    if (secondaryStructureSensitivity) {
                        characterArray[groupIndex] = secStrucState.getType().type;
                    } else {
                        characterArray[groupIndex] = 1;
                    }
                } else {
                    characterArray[groupIndex] = 0;
                }
            }
            groupIndex++;
        }
        return characterArray;
    }


}
