package com.github.pmateusz.bioinformatics.comparing;

import lombok.Getter;

import java.util.ArrayList;

class Node<TValue, TInfo> {

    @Getter
    private final int id;

    @Getter
    private final TValue value;

    @Getter
    private final ArrayList<Edge<TValue, TInfo>> edges = new ArrayList<>();

    Node(int id, TValue value) {
        this.id = id;
        this.value = value;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Node<?, ?> node = (Node<?, ?>) o;

        return id == node.id;

    }

    @Override
    public int hashCode() {
        return id;
    }
}
