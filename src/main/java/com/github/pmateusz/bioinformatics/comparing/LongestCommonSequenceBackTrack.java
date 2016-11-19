package com.github.pmateusz.bioinformatics.comparing;

import javafx.util.Pair;
import lombok.Data;

public class LongestCommonSequenceBackTrack {

    private enum Direction {South, East, SouthEast}

    public static Result Compare(Sequence left, Sequence right, ScoringMatrix scoringMatrix, int indelWeight) {
        Pair<Integer, Direction[][]> result = BackTrack(left, right, scoringMatrix, indelWeight);
        return Output(result, left, right);
    }

    private static Result Output(Pair<Integer, Direction[][]> result, Sequence left, Sequence right) {

        final int score = result.getKey();
        final Direction[][] backtrack = result.getValue();
        int i = backtrack.length - 1;
        int j = backtrack[1].length - 1;

        StringBuilder leftBuilder = new StringBuilder();
        StringBuilder rightBuilder = new StringBuilder();
        while (i >= 0 && j >= 0) {
            switch (backtrack[i][j]) {
                case SouthEast:
                    leftBuilder.append(left.get(i));
                    rightBuilder.append(right.get(j));
                    --i;
                    --j;
                    break;
                case East:
                    leftBuilder.append("-");
                    rightBuilder.append(right.get(j));
                    --j;
                    break;
                case South:
                    leftBuilder.append(left.get(i));
                    rightBuilder.append("-");
                    --i;
                    break;
            }
        }

        while (j >= 0) {
            rightBuilder.append(right.get(j));
            leftBuilder.append("-");
            --j;
        }

        while (i >= 0) {
            rightBuilder.append("-");
            leftBuilder.append(left.get(i));
            --i;
        }

        Sequence leftSequence = new Sequence(leftBuilder.reverse().toString());
        Sequence rightSequence = new Sequence(rightBuilder.reverse().toString());
        return new Result(score, leftSequence, rightSequence);
    }

    private static Pair<Integer, Direction[][]> BackTrack(Sequence left, Sequence right, ScoringMatrix scoringMatrix, int indelWeight) {
        int v = left.length() + 1;
        int w = right.length() + 1;
        int[][] s = new int[v][w];

        int indelWeightToUse = indelWeight * -1;
        Direction[][] backtrack = new Direction[v - 1][w - 1];
        s[0][0] = 0;
        for (int i = 1; i < v; ++i) {
            s[i][0] = s[i - 1][0] + indelWeightToUse;
        }
        for (int j = 1; j < w; ++j) {
            s[0][j] = s[0][j - 1] + indelWeightToUse;
        }

        for (int j = 1; j < w; ++j) {
            for (int i = 1; i < v; ++i) {
                final int goEastScore = s[i][j - 1] + indelWeightToUse;
                final int goSouthScore = s[i - 1][j] + indelWeightToUse;

                final char leftCharacter = left.get(i - 1);
                final char rightCharacter = right.get(j - 1);
                final int goPartSouthEastScore = s[i - 1][j - 1] + scoringMatrix.getWeight(leftCharacter, rightCharacter);

                Direction candidateDirection = Direction.East;
                int candidateScore = goEastScore;

                if (candidateScore <= goSouthScore) {
                    candidateDirection = Direction.South;
                    candidateScore = goSouthScore;
                }

                if (candidateScore <= goPartSouthEastScore) {
                    candidateDirection = Direction.SouthEast;
                    candidateScore = goPartSouthEastScore;
                }

                s[i][j] = candidateScore;
                backtrack[i - 1][j - 1] = candidateDirection;
            }
        }

        return new Pair<>(s[v - 1][w - 1], backtrack);
    }

    @Data
    public static class Result {
        private final int score;
        private final Sequence left;
        private final Sequence right;
    }
}
