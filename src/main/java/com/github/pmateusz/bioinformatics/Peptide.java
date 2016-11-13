package com.github.pmateusz.bioinformatics;

import lombok.Getter;

import java.util.*;
import java.util.stream.Collectors;

public enum Peptide {
    HISTIDINE('H', "His", "Histidine", 137, "CAC", "CAU"),
    GLUTAMINE('Q', "Gin", "Glutamine", 128, "CAA", "CAG"),
    PROLINE('P', "Pro", "Proline", 97, "CCA", "CCC", "CCG", "CCU"),
    ARGININE('R', "Arg", "Arginine", 156, "AGA", "AGG", "CGA", "CGC", "CGG", "CGU"),
    LEUCINE('L', "Leu", "Leucine", 113, "CUA", "CUU", "CUC", "CUG", "UUA", "UUG"),
    ASPARTIC_ACID('D', "Asp", "Aspartic acid", 115, "GAC", "GAU"),
    GLUTAMIC_ACID('E', "Glu", "Glutamic acid", 129, "GAA", "GAG"),
    ALANINE('A', "Ala", "Alanine", 71, "GCA", "GCC", "GCG", "GCU"),
    GLYCINE('G', "Gly", "Glycine", 57, "GGA", "GGC", "GGG", "GGU"),
    VALINE('V', "Val", "Valine", 99, "GUA", "GUC", "GUG", "GUU"),
    TYROSINE('Y', "Tyr", "Tyrosine", 163, "UAC", "UAU"),
    STOP('*', "STP", "STOP", 0, "UAG", "UAA", "UGA"),
    SERINE('S', "Ser", "Serine", 87, "AGC", "AGU", "UCA", "UCC", "UCG", "UCU"),
    CYSTEINE('C', "Cys", "Cysteine", 103, "UGC", "UGU"),
    TRYPTOPHAN('W', "Trp", "Tryptophan", 186, "UGG"),
    PHENYLALANINE('F', "Phe", "Phenylalanine", 147, "UUC", "UUU"),
    ASPARAGINE('N', "Asn", "Asparagine", 114, "AAC", "AAU"),
    LYSINE('K', "Lys", "Lysine", 128, "AAA", "AAG"),
    THREONINE('T', "Thr", "Threonine", 101, "ACA", "ACC", "ACG", "ACU"),
    ISOLEUCINE('I', "Ile", "Isoleucine", 113, "AUA", "AUC", "AUU"),
    METHIONINE('M', "Met", "Methionine", 131, "AUG");

    public static Collection<Peptide> peptides() {
        return Peptides;
    }

    private static final List<Peptide> Peptides = Arrays.stream(Peptide.values()).filter(peptide -> peptide != STOP).collect(Collectors.toList());
    private static final Map<String, Peptide> TripleCodes = new HashMap<>();
    private static final Map<Character, Peptide> Shortcuts = new HashMap<>();

    static {
        for (Peptide peptide : Peptide.values()) {
            for (String code : peptide.codes) {
                TripleCodes.putIfAbsent(code, peptide);
            }
            Shortcuts.put(peptide.shortcut, peptide);
        }
    }

    @Getter
    private char shortcut;

    @Getter
    private String shortName;

    @Getter
    private String name;

    @Getter
    private String[] codes;

    @Getter
    private int mass;

    Peptide(char shortcut, String shortName, String name, int mass, String... codes) {
        this.shortcut = shortcut;
        this.shortName = shortName;
        this.name = name;
        this.mass = mass;
        this.codes = codes;
    }

    public static Peptide get(Character codon) {
        return Shortcuts.get(codon);
    }

    public static Peptide getByCode(String code) {
        return TripleCodes.get(code);
    }

    public static int getMass(String peptideChain) {
        int mass = 0;
        for (char peptide : peptideChain.toCharArray()) {
            mass += Shortcuts.get(peptide).getMass();
        }
        return mass;
    }

    public static Collection<Integer> getSpectrum(String peptideChain) {
        ArrayList<Integer> spectrum = new ArrayList<>(peptideChain.length());
        for (char shortcut : peptideChain.toCharArray()) {
            spectrum.add(Peptide.get(shortcut).getMass());
        }
        return spectrum;
    }

    public static String toNucleotides(String code) {
        final int Length = 3;
        final StringBuilder translationBuilder = new StringBuilder();
        final int length = code.length();
        for (int position = 0; position < length; position += Length) {
            final String codon = code.substring(position, Length);
            translationBuilder.append(Peptide.getByCode(codon));
        }
        return translationBuilder.toString();
    }
}
