package com.github.pmateusz.bioinformatics.comparing;

import javafx.util.Pair;
import lombok.Data;

class GlobalAlignmentProblem {

    static Result compareGlobal(String left, String right, IScoringMatrix scoringMatrix, int indelWeight) {
        Pair<Integer, Direction[][]> result = backTrack(left, right, scoringMatrix, indelWeight);
        return output(result, left, right);
    }

    static GlobalAlignmentProblem.Result overlapAlignmentProblem(String left, String right) {
        final int penalty = -2;
        final int match = 1;
        final int mismatch = -2;

        final IScoringMatrix unityScoringMatrix = (left1, right1) -> {
            if (left1 == right1) {
                return match;
            }
            return mismatch;
        };
        final int rows = left.length() + 1;
        final int columns = right.length() + 1;

        int[][] score = new int[rows][columns];
        Direction[][] movements = new Direction[rows][columns];
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                score[row][column] = 0;
            }
        }

        for (int row = 1; row < rows; ++row) {
            for (int column = 1; column < columns; ++column) {
                int scoreGoEast = score[row][column - 1] + penalty;
                int scoreGoSouth = score[row - 1][column] + penalty;
                int scoreGoSouthEast = score[row - 1][column - 1] + unityScoringMatrix.getScore(left.charAt(row - 1), right.charAt(column - 1));

                Direction candidateMovement = Direction.South;
                int candidateScore = scoreGoSouth;

                if (candidateScore < scoreGoEast) {
                    candidateScore = scoreGoEast;
                    candidateMovement = Direction.East;
                }

                if (candidateScore < scoreGoSouthEast) {
                    candidateScore = scoreGoSouthEast;
                    candidateMovement = Direction.SouthEast;
                }

                score[row][column] = candidateScore;
                movements[row][column] = candidateMovement;
            }
        }

        int candidateColumn = 0;
        final int candidateRow = rows - 1;
        int candidateDistance = score[candidateRow][candidateColumn];
        for (int column = 1; column < columns; ++column) {
            if (candidateDistance < score[candidateRow][column]) {
                candidateDistance = score[candidateRow][column];
                candidateColumn = column;
            }
        }

        int i = candidateRow;
        int j = candidateColumn;
        final StringBuilder leftBuilder = new StringBuilder();
        final StringBuilder rightBuilder = new StringBuilder();
        while (i > 0 && j > 0) {
            switch (movements[i][j]) {
                case East:
                    rightBuilder.append(right.charAt(j - 1));
                    leftBuilder.append('-');
                    --j;
                    break;
                case South:
                    leftBuilder.append(left.charAt(i - 1));
                    rightBuilder.append('-');
                    --i;
                    break;
                case SouthEast:
                    rightBuilder.append(right.charAt(j - 1));
                    leftBuilder.append(left.charAt(i - 1));
                    --i;
                    --j;
                    break;
            }
        }

        return new Result(candidateDistance, leftBuilder.reverse().toString(), rightBuilder.reverse().toString());
    }

    static GlobalAlignmentProblem.Result fittingAlignmentProblem(String left, String right) {
        final int penalty = -1;
        final int match = 1;
        final int mismatch = -1;

        final IScoringMatrix unityScoringMatrix = (left1, right1) -> {
            if (left1 == right1) {
                return match;
            }
            return mismatch;
        };
        final int rows = left.length() + 1;
        final int columns = right.length() + 1;

        int[][] score = new int[rows][columns];
        Direction[][] movements = new Direction[rows][columns];
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                score[row][column] = 0;
            }
        }

        for (int row = 0; row < rows; ++row) {
            score[row][0] = 0;
        }

        for (int column = 0; column < columns; ++column) {
            score[0][column] = column * penalty;
        }

        for (int row = 1; row < rows; ++row) {
            for (int column = 1; column < columns; ++column) {
                int scoreGoEast = score[row][column - 1] + penalty;
                int scoreGoSouth = score[row - 1][column] + penalty;
                int scoreGoSouthEast = score[row - 1][column - 1] + unityScoringMatrix.getScore(left.charAt(row - 1), right.charAt(column - 1));

                Direction candidateMovement = Direction.South;
                int candidateScore = scoreGoSouth;

                if (candidateScore < scoreGoEast) {
                    candidateScore = scoreGoEast;
                    candidateMovement = Direction.East;
                }

                if (candidateScore < scoreGoSouthEast) {
                    candidateScore = scoreGoSouthEast;
                    candidateMovement = Direction.SouthEast;
                }

                score[row][column] = candidateScore;
                movements[row][column] = candidateMovement;
            }
        }

        final int candidateColumn = columns - 1;
        int candidateRow = 0;
        int candidateDistance = score[candidateRow][candidateColumn];
        for (int row = 1; row < rows; ++row) {
            if (candidateDistance < score[row][candidateColumn]) {
                candidateDistance = score[row][candidateColumn];
                candidateRow = row;
            }
        }

        int i = candidateRow;
        int j = candidateColumn;
        final StringBuilder leftBuilder = new StringBuilder();
        final StringBuilder rightBuilder = new StringBuilder();
        while (i > 0 && j > 0) {
            switch (movements[i][j]) {
                case East:
                    rightBuilder.append(right.charAt(j - 1));
                    leftBuilder.append('-');
                    --j;
                    break;
                case South:
                    leftBuilder.append(left.charAt(i - 1));
                    rightBuilder.append('-');
                    --i;
                    break;
                case SouthEast:
                    rightBuilder.append(right.charAt(j - 1));
                    leftBuilder.append(left.charAt(i - 1));
                    --i;
                    --j;
                    break;
            }
        }

        return new Result(candidateDistance, leftBuilder.reverse().toString(), rightBuilder.reverse().toString());
    }

    static int computeEditDistance(String left, String right) {
        final int rows = right.length() + 1;
        final int columns = left.length() + 1;

        int[][] score = new int[rows][columns];
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                score[row][column] = 0;
            }
        }

        for (int row = 0; row < rows; ++row) {
            score[row][0] = row;
        }

        for (int column = 0; column < columns; ++column) {
            score[0][column] = column;
        }

        for (int row = 1; row < rows; ++row) {
            for (int column = 1; column < columns; ++column) {
                final int cost = left.charAt(column - 1) == right.charAt(row - 1) ? 0 : 1;
                int candidateScore = Math.min(score[row - 1][column], score[row][column - 1]) + 1;
                candidateScore = Math.min(candidateScore, score[row - 1][column - 1] + cost);
                score[row][column] = candidateScore;
            }
        }


        return score[rows - 1][columns - 1];
    }

    static Result compareLocal(String left, String right, IScoringMatrix scoringMatrix, int indelWeight) {
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
                final int goPartSouthEastScore = s[i - 1][j - 1] + scoringMatrix.getScore(leftCharacter, rightCharacter);

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

    private static Result output(Pair<Integer, Direction[][]> result, String left, String right) {

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

    private static Pair<Integer, Direction[][]> backTrack(String left, String right, IScoringMatrix scoringMatrix, int indelWeight) {
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
                final int goPartSouthEastScore = s[i - 1][j - 1] + scoringMatrix.getScore(leftCharacter, rightCharacter);

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
