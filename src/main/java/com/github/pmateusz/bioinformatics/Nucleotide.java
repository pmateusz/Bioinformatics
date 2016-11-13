package com.github.pmateusz.bioinformatics;

import lombok.Data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

public enum Nucleotide {
    ADENINE(0, "A", "Adenine"),
    GUANINE(1, "G", "Guanine"),
    CITOSINE(2, "C", "Citosine"),
    THYMINE(3, "T", "Thymine"),
    URACIL(4, "U", "Uracil");

    static {
        ADENINE.complement = CITOSINE;
        CITOSINE.complement = ADENINE;
        GUANINE.complement = THYMINE;
        THYMINE.complement = GUANINE;
        URACIL.complement = GUANINE;
    }

    private final int number;
    private final String shortName;
    private final String longName;
    private Nucleotide complement;

    Nucleotide(int number, String shortName, String longName) {
        this.number = number;
        this.shortName = shortName;
        this.longName = longName;
    }

    public static Nucleotide parseShortName(String shortName) {
        final String shortNameToUse = shortName.toUpperCase(Locale.ROOT);

        switch (shortNameToUse) {
            case "A":
                return ADENINE;
            case "G":
                return GUANINE;
            case "C":
                return CITOSINE;
            case "T":
                return THYMINE;
            case "U":
                return URACIL;
            default:
                throw new IllegalArgumentException(shortName);
        }
    }

    public static int toNumber(String pattern) {
        if (pattern.length() == 0) {
            return 0;
        }

        int result = 0;
        for (char symbol : pattern.toCharArray()) {
            result = 4 * result + toNumber(symbol);
        }
        return result;
    }

    public static int toNumber(char symbol) {
        switch (symbol) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            default:
                throw new IllegalArgumentException(String.valueOf(symbol));
        }
    }

    public static String toSymbol(int number, int k) {
        int currentNumber = number;
        ArrayList<Character> pattern = new ArrayList<>();
        while (k > 0) {
            pattern.add(toSymbol(currentNumber % 4));
            currentNumber = currentNumber / 4;
            --k;
        }
        Collections.reverse(pattern);
        StringBuilder builder = new StringBuilder();
        pattern.forEach(builder::append);
        return builder.toString();
    }

    public static char toSymbol(int number) {
        switch (number) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                throw new IllegalArgumentException(Integer.valueOf(number).toString());
        }
    }
}
