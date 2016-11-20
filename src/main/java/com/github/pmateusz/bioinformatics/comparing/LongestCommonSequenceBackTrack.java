package com.github.pmateusz.bioinformatics.comparing;

import javafx.util.Pair;
import lombok.Data;

class LongestCommonSequenceBackTrack {

    static Result compare(String left, String right, ScoringMatrix scoringMatrix, int indelWeight) {
        Pair<Integer, Direction[][]> result = BackTrack(left, right, scoringMatrix, indelWeight);
        return Output(result, left, right);
    }

    static Result compareLocal(String left, String right, ScoringMatrix scoringMatrix, int indelWeight) {
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

                final char leftCharacter = left.charAt(i - 1);
                final char rightCharacter = right.charAt(j - 1);
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

                if (candidateScore < 0) {
                    candidateDirection = null;
                    candidateScore = 0;
                }

                s[i][j] = candidateScore;
                backtrack[i - 1][j - 1] = candidateDirection;
            }
        }

        int candidateScore = Integer.MIN_VALUE;
        int candidateRow = 0, candidateColumn = 0;
        for (int row = 1; row < v; ++row) {
            for (int column = 1; column < w; ++column) {
                if (candidateScore < s[row][column]) {
                    candidateScore = s[row][column];
                    candidateRow = row;
                    candidateColumn = column;
                }
            }
        }

        StringBuilder leftBuilder = new StringBuilder();
        StringBuilder rightBuilder = new StringBuilder();
        int row = candidateRow - 1;
        int column = candidateColumn - 1;
        while (row >= 0 && column >= 0) {
            Direction candidateMove = backtrack[row][column];
            if (candidateMove == null) {
                break;
            }
            switch (backtrack[row][column]) {
                case SouthEast:
                    leftBuilder.append(left.charAt(row));
                    rightBuilder.append(right.charAt(column));
                    --row;
                    --column;
                    break;
                case East:
                    leftBuilder.append("-");
                    rightBuilder.append(right.charAt(column));
                    --column;
                    break;
                case South:
                    leftBuilder.append(left.charAt(row));
                    rightBuilder.append("-");
                    --row;
                    break;
            }
        }

        String leftResult = leftBuilder.reverse().toString();
        String rightResult = rightBuilder.reverse().toString();
        return new Result(candidateScore, leftResult, rightResult);
    }

    private static Result Output(Pair<Integer, Direction[][]> result, String left, String right) {

        final int score = result.getKey();
        final Direction[][] backtrack = result.getValue();
        int i = backtrack.length - 1;
        int j = backtrack[1].length - 1;

        StringBuilder leftBuilder = new StringBuilder();
        StringBuilder rightBuilder = new StringBuilder();
        while (i >= 0 && j >= 0) {
            switch (backtrack[i][j]) {
                case SouthEast:
                    leftBuilder.append(left.charAt(i));
                    rightBuilder.append(right.charAt(j));
                    --i;
                    --j;
                    break;
                case East:
                    leftBuilder.append("-");
                    rightBuilder.append(right.charAt(j));
                    --j;
                    break;
                case South:
                    leftBuilder.append(left.charAt(i));
                    rightBuilder.append("-");
                    --i;
                    break;
            }
        }

        while (j >= 0) {
            rightBuilder.append(right.charAt(j));
            leftBuilder.append("-");
            --j;
        }

        while (i >= 0) {
            rightBuilder.append("-");
            leftBuilder.append(left.charAt(i));
            --i;
        }

        String leftSequence = leftBuilder.reverse().toString();
        String rightSequence = rightBuilder.reverse().toString();
        return new Result(score, leftSequence, rightSequence);
    }

    private static Pair<Integer, Direction[][]> BackTrack(String left, String right, ScoringMatrix scoringMatrix, int indelWeight) {
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

                final char leftCharacter = left.charAt(i - 1);
                final char rightCharacter = right.charAt(j - 1);
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
    static class Result {
        private final int score;
        private final String left;
        private final String right;
    }
}
