package com.github.pmateusz.bioinformatics.comparing;

import lombok.Getter;
import lombok.Setter;

class Edge<TValue, TInfo> {

    @Getter
    private final int id;

    @Getter
    private final int weight;

    @Getter
    private final TInfo info;

    @Getter
    @Setter
    private Node<TValue, TInfo> source;

    @Getter
    @Setter
    private Node<TValue, TInfo> destination;

    Edge(int id, int weight, TInfo info) {
        this.id = id;
        this.weight = weight;
        this.info = info;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Edge<?, ?> edge = (Edge<?, ?>) o;
        return id == edge.id;

    }

    @Override
    public int hashCode() {
        return id;
    }

    @Override
    public String toString() {
        String label = "Edge(" + id + ", " + weight + ")";
        if (info != null) {
            label += ' ' + info.toString();
        }
        return label;
    }
}
