package com.github.pmateusz.bioinformatics.comparing;

import javafx.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

public class Graph<TValue, TInfo> {

    private final EdgeFactory edgeFactory = new EdgeFactory();
    private final NodeFactory nodeFactory = new NodeFactory();

    private final ArrayList<Node<TValue, TInfo>> nodes = new ArrayList<>();
    private final ArrayList<Edge<TValue, TInfo>> edges = new ArrayList<>();

    public Graph() {
    }

    public Node<TValue, TInfo> addNode(TValue value) {
        Node<TValue, TInfo> node = nodeFactory.createNode(value);
        nodes.add(node);
        return node;
    }

    public Optional<Node<TValue, TInfo>> findNode(TValue value) {
        return nodes.stream().filter(n -> n.getValue().equals(value)).findFirst();
    }

    public Optional<Node<TValue, TInfo>> findNode(int id) {
        return nodes.stream().filter(n -> n.getId() == id).findFirst();
    }

    public Edge<TValue, TInfo> addEdge(Node<TValue, TInfo> source, Node<TValue, TInfo> destination, int weight, TInfo info) {
        Edge<TValue, TInfo> edge = edgeFactory.createEdge(weight, info);
        edge.setSource(source);
        edge.setDestination(destination);
        edges.add(edge);
        return edge;
    }

    public Pair<Map<Node<TValue, TInfo>, Edge<TValue, TInfo>>, Map<Node<TValue, TInfo>, Integer>> bellmanFord(Node<TValue, TInfo> source) {
        Map<Node<TValue, TInfo>, Integer> distances = new HashMap<>();
        Map<Node<TValue, TInfo>, Edge<TValue, TInfo>> predecessors = new HashMap<>();

        for (Node<TValue, TInfo> node : nodes) {
            distances.put(node, Integer.MAX_VALUE);
            predecessors.put(node, null);
        }

        distances.put(source, 0);

        final int nodesCount = nodes.size();
        for (int i = 0; i < nodesCount; ++i) {
            boolean relaxed = false;
            for (Edge<TValue, TInfo> edge : edges) {
                final Node<TValue, TInfo> edgeSource = edge.getSource();
                final Node<TValue, TInfo> edgeDestination = edge.getDestination();

                if (distances.get(edgeSource) == Integer.MAX_VALUE) {
                    continue;
                }

                final Integer candidateDistance = distances.get(edgeSource) + edge.getWeight();
                if (candidateDistance < distances.get(edgeDestination)) {
                    distances.put(edgeDestination, candidateDistance);
                    predecessors.put(edgeDestination, edge);
                    relaxed = true;
                }
            }

            if (!relaxed) {
                break;
            }
        }

        for (int i = 0; i < nodesCount; ++i) {
            for (Edge<TValue, TInfo> edge : edges) {
                final Node<TValue, TInfo> edgeSource = edge.getSource();
                final Node<TValue, TInfo> edgeDestination = edge.getDestination();
                final Integer candidateDistance = distances.get(edgeSource) + edge.getWeight();
                if (candidateDistance < distances.get(edgeDestination)) {
                    throw new IllegalStateException("Graph contains cycles");
                }
            }
        }

        return new Pair<>(predecessors, distances);
    }

    private class EdgeFactory {

        private int currentId;

        EdgeFactory() {
            this.currentId = 0;
        }

        public Edge<TValue, TInfo> createEdge(int weight, TInfo info) {
            return new Edge<>(currentId++, weight, info);
        }
    }

    private class NodeFactory {

        private int currentId;

        NodeFactory() {
            this.currentId = 0;
        }

        public Node<TValue, TInfo> createNode(TValue value) {
            return new Node<>(currentId++, value);
        }
    }
}
