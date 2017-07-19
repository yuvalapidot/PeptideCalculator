package model;

/**
 * Created by yuvalapidot.
 */
public class ScoreBoard {

    public int[] primaryScore;
    public int[] secondaryScore;
    public int[] peptideSize;
    public int[] link;

    public ScoreBoard(int chainSize, int peptideSize) {
        this.primaryScore = new int[chainSize];
        this.secondaryScore = new int[chainSize];
        this.peptideSize = new int[chainSize];
        this.link = new int[chainSize];
        initialize(chainSize, peptideSize);
    }

    private void initialize(int chainSize, int peptideSize) {
        for (int i = 1; i <= peptideSize; i++) {
            this.primaryScore[chainSize - i] = 0;
            this.secondaryScore[chainSize - i] = peptideSize - i;
            this.peptideSize[chainSize - i] = i;
            this.link[chainSize - i] = -1;
        }
    }
}
