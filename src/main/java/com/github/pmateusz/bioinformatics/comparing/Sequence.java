package com.github.pmateusz.bioinformatics.comparing;

public class Sequence {

    private final String sequence;
    private final char[] rawSequence;

    public Sequence(String sequence) {
        this.sequence = sequence;
        this.rawSequence = sequence.toCharArray();
    }

    public int length() {
        return sequence.length();
    }

    public char get(int index) {
        return rawSequence[index];
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final Sequence other = (Sequence) o;
        return !(sequence != null ? !sequence.equals(other.sequence) : other.sequence != null);

    }

    @Override
    public int hashCode() {
        return sequence != null ? sequence.hashCode() : 0;
    }

    @Override
    public String toString() {
        return sequence;
    }
}
