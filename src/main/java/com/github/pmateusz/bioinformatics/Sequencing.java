package com.github.pmateusz.bioinformatics;

import com.google.common.collect.Iterables;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Sequencing {

    static Collection<Integer> spectralConvolution(Collection<Integer> spectrum) {
        ArrayList<Integer> spectrumToUse = new ArrayList<>(spectrum);
        Collections.sort(spectrumToUse);

        if (spectrumToUse.size() == 0) {
            return Collections.emptyList();
        }

        if (spectrumToUse.get(0) != 0) {
            ArrayList<Integer> localList = new ArrayList<>(spectrumToUse.size() + 1);
            localList.add(0);
            localList.addAll(spectrumToUse);
            spectrumToUse = localList;
        }

        Map<Integer, Integer[]> convolution = new HashMap<>();
        for (Integer key : spectrumToUse) {
            Stream<Integer> spectrumStream = spectrum.stream()
                    .filter(m -> m < key)
                    .map(m -> key - m);
            if (convolution.containsKey(key)) {
                spectrumStream = Stream.concat(spectrumStream, Arrays.stream(convolution.get(key)));
            }
            convolution.put(key, spectrumStream.toArray(Integer[]::new));
        }

        Map<Integer, Integer> groupByMass = new HashMap<>();
        for (Integer[] masses : convolution.values()) {
            for (Integer mass : masses) {
                Integer frequency = groupByMass.getOrDefault(mass, 0) + 1;
                groupByMass.put(mass, frequency);
            }
        }

        final int CutOffThreshold = 0;
        return groupByMass.entrySet().stream()
                .filter(pair -> pair.getValue() > CutOffThreshold && pair.getKey() > 0)
                .map(pair -> Collections.nCopies(pair.getValue(), pair.getKey()))
                .flatMap(Collection::stream)
                .sorted()
                .collect(Collectors.toList());
    }

    public static PeptideChain convolutionCyclopeptideSeqeuencing(Collection<Integer> spectrum, int M, int N) {
        final int MIN_PEPTIDE_MASS = 57;
        final int MAX_PEPTIDE_MASS = 200;

        Collection<Integer> convolutedSpectrum = spectralConvolution(spectrum);
        Map<Integer, Integer> spectrumFrequency = new HashMap<>();
        for (Integer mass : convolutedSpectrum) {
            if (mass >= MIN_PEPTIDE_MASS && mass <= MAX_PEPTIDE_MASS) {
                final Integer frequency = spectrumFrequency.getOrDefault(mass, 0) + 1;
                spectrumFrequency.put(mass, frequency);
            }
        }

        LinkedHashMap<Integer, Integer> sortedFrequencySpectrum = spectrumFrequency.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByValue(Collections.reverseOrder()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (k, v) -> v, LinkedHashMap::new));

        final int cutOffToUse = Iterables.get(sortedFrequencySpectrum.values(), Math.min(sortedFrequencySpectrum.size() - 1, M));
        long massesToExclude = sortedFrequencySpectrum.entrySet().stream().filter(pair -> pair.getValue() < cutOffToUse).count();
        List<Integer> peptides = Stream.concat(sortedFrequencySpectrum.entrySet()
                .stream()
                .limit(sortedFrequencySpectrum.size() - massesToExclude)
                .map(Map.Entry::getKey), Peptide.peptides().stream().map(Peptide::getMass))
                .distinct()
                .sorted()
                .collect(Collectors.toList());

        List<Integer> spectrumToUse = spectrum.stream().sorted().collect(Collectors.toList());
        return leaderBoardCyclopeptideSequencing(spectrumToUse, peptides, M);
    }

    private static PeptideChain leaderBoardCyclopeptideSequencing(Collection<Integer> spectrum, Collection<Integer> peptides, int N) {
        PeptideChain leaderPeptide = PeptideChain.empty();
        HashSet<PeptideChain> leaderBoard = new HashSet<>();
        leaderBoard.add(leaderPeptide);

        int parentMass = spectrum.stream().max(Comparator.naturalOrder()).orElse(0);
        while (!leaderBoard.isEmpty()) {
            HashSet<PeptideChain> nextLeaderBoard = new HashSet<>();
            for (PeptideChain peptideChain : leaderBoard) {
                for (Integer peptide : peptides) {
                    nextLeaderBoard.add(peptideChain.add(peptide));
                }
            }
            leaderBoard = nextLeaderBoard;

            Iterator<PeptideChain> leaderBoardIterator = leaderBoard.iterator();
            while (leaderBoardIterator.hasNext()) {
                final PeptideChain peptideChain = leaderBoardIterator.next();
                int peptideMass = peptideChain.getTotalMass();
                if (peptideMass == parentMass) {
                    if (peptideChain.getLinearScore(spectrum) > leaderPeptide.getLinearScore(spectrum)) {
                        leaderPeptide = peptideChain;
                    }
                } else if (peptideMass > parentMass) {
                    leaderBoardIterator.remove();
                }
            }

            leaderBoard = new HashSet<>(trimLeaderBoard(leaderBoard, spectrum, N));
        }

        return leaderPeptide;
    }

    private static Collection<PeptideChain> trimLeaderBoard(HashSet<PeptideChain> leaderBoard, Collection<Integer> spectrum, int N) {
        HashMap<PeptideChain, Integer> scores = new HashMap<>();
        for (PeptideChain peptide : leaderBoard) {
            scores.put(peptide, peptide.getScore(spectrum));
        }

        List<PeptideChain> peptidesInOrder = leaderBoard.stream().sorted((left, right) -> scores.get(right) - scores.get(left)).collect(Collectors.toList());
        int leaderBoardSize = leaderBoard.size();
        if (leaderBoardSize <= N) {
            return leaderBoard;
        }

        int cutoff = scores.get(peptidesInOrder.get(N));
        for (int j = N + 1; j < leaderBoardSize; ++j) {
            final PeptideChain peptideChain = peptidesInOrder.get(j);
            if (scores.get(peptideChain) < cutoff) {
                return peptidesInOrder.stream().limit(j).collect(Collectors.toList());
            }
        }
        return peptidesInOrder;
    }

    public static Collection<PeptideChain> cyclopeptideSequencing(Collection<Integer> spectrum) {
        HashSet<PeptideChain> results = new HashSet<>();
        HashSet<Integer> spectrumPeptides = new HashSet<>(spectrum);
        PeptideChain spectrumChain = new PeptideChain(new ArrayList<>(spectrum));

        int parentMass = spectrum.stream().max(Comparator.naturalOrder()).orElse(0);
        HashSet<PeptideChain> peptideChains = new HashSet<>();
        peptideChains.add(PeptideChain.empty());

        while (!peptideChains.isEmpty()) {
            HashSet<PeptideChain> nextPeptideChains = new HashSet<>();
            for (PeptideChain peptide : peptideChains) {
                for (Peptide component : Peptide.peptides()) {
                    nextPeptideChains.add(peptide.add(component.getMass()));
                }
            }
            peptideChains = nextPeptideChains;

            Iterator<PeptideChain> peptideIterator = peptideChains.iterator();
            while (peptideIterator.hasNext()) {
                final PeptideChain peptideChain = peptideIterator.next();
                if (peptideChain.getTotalMass() == parentMass) {
                    PeptideChain testSpectrum = peptideChain.getCyclicSpectrum();
                    if (testSpectrum.equals(spectrumChain)) {
                        results.add(peptideChain);
                    }

                    peptideIterator.remove();
                } else if (!spectrumPeptides.contains(peptideChain.getTotalMass())) {
                    peptideIterator.remove();
                }
            }
        }

        return results;
    }
}
