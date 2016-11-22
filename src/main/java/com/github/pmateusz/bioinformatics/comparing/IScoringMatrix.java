package com.github.pmateusz.bioinformatics.comparing;

public interface IScoringMatrix {

    int getScore(Character left, Character right);
}
