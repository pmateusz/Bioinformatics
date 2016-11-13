package com.github.pmateusz.bioinformatics;

import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;
import lombok.Getter;

import java.lang.reflect.Array;
import java.util.*;
import java.util.stream.Collectors;

public class PeptideChain implements Iterable<Integer> {
    private final ArrayList<Integer> massChain;
    @Getter
    private final int totalMass;

    public static PeptideChain empty() {
        return new PeptideChain(new ArrayList<>());
    }

    public static PeptideChain fromCode(String peptideChain) {
        ArrayList<Integer> massChain = new ArrayList<>(peptideChain.length());
        for (char character : peptideChain.toCharArray()) {
            final Integer mass = Peptide.get(character).getMass();
            massChain.add(mass);
        }

        return new PeptideChain(massChain);
    }

    public PeptideChain(ArrayList<Integer> massChain) {
        this.massChain = massChain;
        this.totalMass = massChain.stream().collect(Collectors.summingInt(p -> p));
    }

    public PeptideChain add(Integer mass) {
        ArrayList<Integer> newMassChain = new ArrayList<>(massChain.size() + 1);
        newMassChain.addAll(massChain);
        newMassChain.add(mass);
        return new PeptideChain(newMassChain);
    }

    public Collection<Integer> getMassChain() {
        return Collections.unmodifiableCollection(massChain);
    }

    public PeptideChain getLinearSpectrum() {
        ArrayList<Integer> prefixMass = new ArrayList<>();
        prefixMass.add(0);

        int lastPrefixMass = 0;
        for(Integer mass: massChain) {
            lastPrefixMass += mass;
            prefixMass.add(lastPrefixMass);
        }

        ArrayList<Integer> linearSpectrum = new ArrayList<>();
        linearSpectrum.add(0);

        final int peptideLength = length();
        for (int i = 0; i < peptideLength; ++i) {
            for (int j = i + 1; j <= peptideLength; ++j) {
                int mass = prefixMass.get(j) - prefixMass.get(i);
                linearSpectrum.add(mass);
            }
        }

        Collections.sort(linearSpectrum);
        return new PeptideChain(linearSpectrum);
    }

    public PeptideChain getCyclicSpectrum() {
        ArrayList<Integer> prefixMass = new ArrayList<>();
        prefixMass.add(0);

        int lastPrefixMass = 0;
        for (Integer mass : massChain) {
            lastPrefixMass += mass;
            prefixMass.add(lastPrefixMass);
        }

        final int peptideMass = Iterables.getLast(prefixMass);
        ArrayList<Integer> cyclicSpectrum = new ArrayList<>();
        cyclicSpectrum.add(0);

        final int peptideLength = length();
        for (int i = 0; i < peptideLength; ++i) {
            for (int j = i + 1; j <= peptideLength; ++j) {
                int mass = prefixMass.get(j) - prefixMass.get(i);
                cyclicSpectrum.add(mass);
                if (i > 0 && j < peptideLength) {
                    mass = peptideMass - mass;
                    cyclicSpectrum.add(mass);
                }
            }
        }

        Collections.sort(cyclicSpectrum);
        return new PeptideChain(cyclicSpectrum);
    }

    public int getScore(Collection<Integer> spectrum) {
        PeptideChain cyclicPeptideChain = getCyclicSpectrum();
        return getScore(cyclicPeptideChain, spectrum);
    }

    public int getLinearScore(Collection<Integer> spectrum) {
        PeptideChain linearPeptideChain = getLinearSpectrum();
        return getScore(linearPeptideChain, spectrum);
    }

    private static int getScore(PeptideChain peptideChain, Collection<Integer> spectrum) {
        HashMap<Integer, Integer> frequency = new HashMap<>();
        for (Integer mass : peptideChain) {
            int count = frequency.getOrDefault(mass, 0);
            frequency.put(mass, count + 1);
        }

        int score = 0;
        for (Integer mass : spectrum) {
            if (frequency.containsKey(mass)) {
                int count = frequency.get(mass);
                if (count > 0) {
                    ++score;
                    frequency.put(mass, --count);
                }
            }
        }
        return score;
    }

    public int length() {
        return massChain.size();
    }

    public String toString() {
        return Joiner.on("-").join(massChain);
    }

    public int hashCode() {
        return massChain.hashCode();
    }

    public boolean equals(Object other) {
        return this == other
                || other instanceof PeptideChain
                && massChain.equals(((PeptideChain) other).massChain);

    }

    @Override
    public Iterator<Integer> iterator() {
        return massChain.iterator();
    }
}
