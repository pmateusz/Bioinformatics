package com.github.pmateusz.bioinformatics.comparing;

import com.google.common.base.Charsets;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class ScoringMatrix {

    private final Map<Character, Integer> indexMap;
    private final int[][] weightMatrix;

    public ScoringMatrix(Map<Character, Integer> indexMap, int[][] weightMatrix) {
        this.indexMap = indexMap;
        this.weightMatrix = weightMatrix;
    }

    public int getWeight(Character left, Character right) {
        final int row = indexMap.get(left);
        final int column = indexMap.get(right);
        return weightMatrix[row][column];
    }

    public static ScoringMatrix load(String resourcePath) throws IOException {
        ClassLoader classLoader = ScoringMatrix.class.getClassLoader();
        URL resourceUrl = classLoader.getResource(resourcePath);
        if (resourceUrl == null) {
            throw new IllegalArgumentException(resourcePath);
        }

        String rawFilePath = resourceUrl.getFile();
        // Paths.get() fails to parse path if it starts with '/', for example /C:/workspace
        if (rawFilePath.startsWith("/")) {
            rawFilePath = rawFilePath.substring(1);
        }

        Path pathToUse = Paths.get(rawFilePath);
        List<String> lines = Files.readAllLines(pathToUse, Charsets.UTF_8);

        List<Character> characters = Arrays.stream(lines.get(0).split("\\s+"))
                .filter(pattern -> pattern.length() > 0)
                .map(singleLetter -> singleLetter.charAt(0))
                .collect(Collectors.toList());
        final HashMap<Character, Integer> indexMap = new HashMap<>();
        int index = 0;
        for (Character character : characters) {
            indexMap.put(character, index++);
        }

        final int length = characters.size();
        final int[][] scoreMatrix = new int[length][length];
        for (int row = 0; row < length; ++row) {
            List<Integer> columnWeights = Arrays.stream(lines.get(row + 1).split("\\s+"))
                    .skip(1)
                    .map(Integer::parseInt)
                    .collect(Collectors.toList());
            int column = 0;
            for (Integer weight : columnWeights) {
                scoreMatrix[row][column++] = weight;
            }
        }

        return new ScoringMatrix(indexMap, scoreMatrix);
    }
}
