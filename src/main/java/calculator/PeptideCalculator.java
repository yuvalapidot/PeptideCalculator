package calculator;

import model.Peptide;
import model.ScoreBoard;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.secstruc.SecStrucState;
import org.biojava.nbio.structure.secstruc.SecStrucType;

import java.util.*;

public class PeptideCalculator {

    private static final String SECSTRUCT_PROPERTY_KEY = "secstruc";

    public List<Map<String, List<Peptide>>> structureToPeptides(Structure structure, int peptideBaseSize,
                                                              int peptideMinSize, int peptideMaxSize,
                                                              int peptideOverlapBase, int peptideOverlapMax,
                                                              Set<SecStrucType> secondaryStructureTypes,
                                                              boolean secondaryStructureSensitivity) throws StructureException {
        assignSecStructureProperty(structure);
        List<Map<String, List<Peptide>>> models = new ArrayList<Map<String, List<Peptide>>>();
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

    private Map<String, List<Peptide>> modelToPeptides(List<Chain> model, int peptideBaseSize, int peptideMinSize, int peptideMaxSize, int peptideOverlapBase, int peptideOverlapMax, Set<SecStrucType> secondaryStructureTypes, boolean secondaryStructureSensitivity) {
        Map<String, List<Peptide>> chainMap = new HashMap<String, List<Peptide>>();
        for (Chain chain : model) {
            List<Peptide> peptides = chainToPeptides(chain, peptideBaseSize, peptideMinSize,
                    peptideMaxSize, peptideOverlapBase, peptideOverlapMax,
                    secondaryStructureTypes, secondaryStructureSensitivity);
            chainMap.put(chain.getChainID(), peptides);
        }
        return chainMap;
    }

    private List<Peptide> chainToPeptides(Chain chain, int peptideBaseSize, int peptideMinSize, int peptideMaxSize, int peptideOverlapBase, int peptideOverlapMax, Set<SecStrucType> secondaryStructureTypes, boolean secondaryStructureSensitivity) {
        char[] characterArray = formatChainAsCharacterArray(chain, secondaryStructureTypes, secondaryStructureSensitivity);
        ScoreBoard scoreBoard = calculateScoreBoard(characterArray, peptideBaseSize, peptideMinSize, peptideMaxSize, peptideOverlapBase, peptideOverlapMax);
        List<Integer> peptideIndexChain = calculatePeptideIndexChain(scoreBoard);
        return createPeptides(scoreBoard, chain, peptideIndexChain);
    }

    private List<Peptide> createPeptides(ScoreBoard scoreBoard, Chain chain, List<Integer> peptideIndexChain) {
        List<Peptide> peptides = new ArrayList<Peptide>();
        for (int index : peptideIndexChain) {
            List<Group> aminoAcids = chain.getSeqResGroups().subList(index, index + scoreBoard.peptideSize[index]);
            Peptide peptide = new Peptide(aminoAcids, index, scoreBoard.peptideSize[index]);
            peptides.add(peptide);
        }
        return peptides;
    }

    private List<Integer> calculatePeptideIndexChain(ScoreBoard scoreBoard) {
        List<Integer> peptideIndexChain = new ArrayList<Integer>();
        int next = 0;
        while (next >= 0) {
            peptideIndexChain.add(next);
            next = scoreBoard.link[next];
        }
        return peptideIndexChain;
    }

    private ScoreBoard calculateScoreBoard(char[] characterArray, int peptideBaseSize, int peptideMinSize, int peptideMaxSize, int peptideOverlapBase, int peptideOverlapMax) {
        ScoreBoard scoreBoard = new ScoreBoard(characterArray.length, peptideBaseSize);
        int[] sizes = getPossibleSizes(peptideBaseSize, peptideMinSize, peptideMaxSize);
        int[] overlaps = getPossibleOverlaps(peptideOverlapBase, peptideOverlapMax);
        for (int i = characterArray.length - peptideBaseSize - 1; i >= 0; i--) {
            calculateScoreboardColumn(scoreBoard, i, characterArray, sizes, overlaps, peptideBaseSize, peptideOverlapBase);
        }
        return scoreBoard;
    }

    private void calculateScoreboardColumn(ScoreBoard scoreBoard, int index, char[] characterArray, int[] sizes, int[] overlaps, int peptideBaseSize, int peptideOverlapBase) {
        int totalSize = characterArray.length;
        scoreBoard.primaryScore[index] = totalSize;
        scoreBoard.secondaryScore[index] = totalSize;
        for (int size : sizes) {
            for (int overlap : overlaps) {
                int nextPeptide = nextColumn(index, size, overlap, characterArray.length);
                int primaryScore;
                if (nextPeptide < -1) {
                    continue;
                }
                if (nextPeptide == -1) {
                    primaryScore = 0;
                } else {
                    primaryScore = calculatePrimaryScore(characterArray, index, nextPeptide, size, scoreBoard);
                }
                if (primaryScore > scoreBoard.primaryScore[index]) {
                    continue;
                }
                int secondaryScore = calculateSecondaryScore(size, overlap, peptideBaseSize, peptideOverlapBase);
                if ((primaryScore < scoreBoard.primaryScore[index]) || (secondaryScore < scoreBoard.secondaryScore[index])) {
                    replaceBestPeptide(index, nextPeptide, primaryScore, secondaryScore, size, scoreBoard);
                }
            }
        }
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

    private int[] getPossibleOverlaps(int peptideOverlapBase, int peptideOverlapMax) {
        int min = peptideOverlapBase;
        int max = Math.max(peptideOverlapBase, peptideOverlapMax);
        int[] overlaps = new int[max - min + 1];
        int index = 0;
        for (int i = min; i <= max; i++) {
            overlaps[index++] = i;
        }
        return overlaps;
    }

    private int[] getPossibleSizes(int peptideBaseSize, int peptideMinSize, int peptideMaxSize) {
        int min = Math.min(peptideBaseSize, peptideMinSize);
        int max = Math.max(peptideBaseSize, peptideMaxSize);
        int[] sizes = new int[max - min + 1];
        int index = 0;
        for (int i = min; i <= max; i++) {
            sizes[index++] = i;
        }
        return sizes;
    }

    private int nextColumn(int column, int size, int overlap, int numberOfColumns) {
        if (column + size > numberOfColumns) {
            return -2;
        }
        if (overlap >= size) {
            return -2;
        }
        if (column + size == numberOfColumns) {
            return -1;
        }
        return column + size - overlap;
    }

    private int calculatePrimaryScore(char[] characterArray, int column, int nextPeptide, int size, ScoreBoard scoreBoard) {
        int primaryScore = scoreBoard.primaryScore[nextPeptide];
        if (isCrossingPeptide(characterArray, column, size, nextPeptide)) {
            primaryScore++;
        }
        return primaryScore;
    }

    private boolean isCrossingPeptide(char[] characterArray, int column, int size, int nextPeptide) {
        char lastIn = characterArray[column + size - 1];
        char firstOut = characterArray[column + size];
        if ((lastIn > 0) && (lastIn == firstOut)) {
            int sequenceStart = calculateSecondarySequenceStart(characterArray, column, lastIn);
            if (nextPeptide <= sequenceStart) {
                return false;
            }
            return true;
        }
        return false;
    }

    private int calculateSecondarySequenceStart(char[] characterArray, int column, char lastIn) {
        while (characterArray[column] == lastIn) {
            column--;
        }
        return column + 1;
    }

    private int calculateSecondaryScore(int size, int overlap, int peptideBaseSize, int peptideOverlapBase) {
        return Math.abs(peptideBaseSize - size) + Math.abs(peptideOverlapBase - overlap);
    }

    private void replaceBestPeptide(int index, int nextPeptide, int primaryScore, int secondaryScore, int size, ScoreBoard scoreBoard) {
        scoreBoard.primaryScore[index] = primaryScore;
        scoreBoard.secondaryScore[index] = secondaryScore;
        scoreBoard.peptideSize[index] = size;
        scoreBoard.link[index] = nextPeptide;
    }
}
