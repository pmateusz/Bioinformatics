package com.github.pmateusz.bioinformatics.comparing;

import javafx.util.Pair;
import lombok.Data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

public class LocalAlignmentProblem {

    @Data
    private static class DirectionInfo {

        private final Character leftLabel;
        private final Character rightLabel;
        private final Direction direction;

        DirectionInfo(Character leftLabel, Character rightLabel, Direction direction) {
            this.leftLabel = leftLabel;
            this.rightLabel = rightLabel;
            this.direction = direction;
        }

        @Override
        public String toString() {
            return "(" + rightLabel + ", " + leftLabel + ")";
        }
    }

    @Data
    public static class Result {
        private final int score;
        private final String left;
        private final String right;
    }

    public static Result compare(String left, String right, ScoringMatrix scoringMatrix, int indelWeight) {
        final Graph<Void, DirectionInfo> graph = new Graph<>();

        final int rows = left.length() + 1;
        final int columns = right.length() + 1;
        final Node<Void, DirectionInfo>[][] nodes = new Node[rows][columns];
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                nodes[row][column] = graph.addNode(null);
            }
        }

        final Character Blank = Character.MIN_VALUE;

        for (int row = 0; row < rows; ++row) {
            for (int column = 1; column < columns; ++column) {
                final Character rightLabel = right.charAt(column - 1);
                final Node<Void, DirectionInfo> source = nodes[row][column - 1];
                final Node<Void, DirectionInfo> destination = nodes[row][column];
                graph.addEdge(source, destination, indelWeight, new DirectionInfo(Blank, rightLabel, Direction.East));
            }
        }

        for (int row = 1; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                final Character leftLabel = left.charAt(row - 1);
                final Node<Void, DirectionInfo> source = nodes[row - 1][column];
                final Node<Void, DirectionInfo> destination = nodes[row][column];
                graph.addEdge(source, destination, indelWeight, new DirectionInfo(leftLabel, Blank, Direction.South));
            }
        }

        for (int row = 1; row < rows; ++row) {
            for (int column = 1; column < columns; ++column) {
                final Character leftLabel = left.charAt(row - 1);
                final Character rightLabel = right.charAt(column - 1);
                final Node<Void, DirectionInfo> source = nodes[row - 1][column - 1];
                final Node<Void, DirectionInfo> destination = nodes[row][column];
                final int weight = scoringMatrix.getWeight(leftLabel, rightLabel) * -1;
                graph.addEdge(source, destination, weight, new DirectionInfo(leftLabel, rightLabel, Direction.SouthEast));
            }
        }

        final Node<Void, DirectionInfo> source = nodes[0][0];
        final Node<Void, DirectionInfo> destination = nodes[rows - 1][columns - 1];
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                final Node<Void, DirectionInfo> currentNode = nodes[row][column];
                if (currentNode != source) {
                    graph.addEdge(source, currentNode, 0, null);
                }
                if (currentNode != destination && currentNode != source) {
                    graph.addEdge(currentNode, destination, 0, null);
                }
            }
        }

        final Pair<Map<Node<Void, DirectionInfo>, Edge<Void, DirectionInfo>>, Map<Node<Void, DirectionInfo>, Integer>> result = graph.bellmanFord(source);
        final Map<Node<Void, DirectionInfo>, Edge<Void, DirectionInfo>> predecessors = result.getKey();
        final Map<Node<Void, DirectionInfo>, Integer> distances = result.getValue();
        final ArrayList<DirectionInfo> path = new ArrayList<>();
        Edge<Void, DirectionInfo> currentEdge = predecessors.get(destination);
        while (currentEdge != null) {
            final DirectionInfo info = currentEdge.getInfo();
            if (info != null) {
                path.add(info);
            }

            currentEdge = predecessors.get(currentEdge.getSource());
        }

        Collections.reverse(path);
        final StringBuilder leftBuilder = new StringBuilder();
        final StringBuilder rightBuilder = new StringBuilder();
        for (DirectionInfo info : path) {
            final Character leftLabel = info.getLeftLabel();
            final Character leftLabelToUse = leftLabel == Blank ? '-' : leftLabel;
            leftBuilder.append(leftLabelToUse);

            final Character rightLabel = info.getRightLabel();
            final Character rightLabelToUse = rightLabel == Blank ? '-' : rightLabel;
            rightBuilder.append(rightLabelToUse);
        }

        return new Result(distances.get(destination) * -1, leftBuilder.toString(), rightBuilder.toString());
    }
}
