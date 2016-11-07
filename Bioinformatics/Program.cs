using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Bioinformatics
{
    class Program
    {
        static void Main (string[] args)
        {
            var lines = File.ReadAllLines (@"C:\Users\matt\Downloads\dataset_102_8.txt");
            var N = int.Parse (lines.First ());
            var input = lines.Skip (1).First ().Split (' ').Select (int.Parse).ToArray ();
            var result = LeaderboardCyclopeptideSequencing (input, N);
            result = string.Join ("-", result.Select (p => Patterns.Peptides[p].ToString ()));
            var formattedResult = new[] {
                string.Join(Environment.NewLine, result)
            };
            File.WriteAllLines (@"C:\Users\matt\Desktop\output.txt", formattedResult);
            Console.WriteLine (string.Join (Environment.NewLine, formattedResult));
            Console.ReadLine ();
        }

        private static string LeaderboardCyclopeptideSequencing (IEnumerable<int> spectrum, int N)
        {
            var leaderboard = new HashSet<string> ();
            leaderboard.Add (string.Empty);

            var parentMass = spectrum.Max ();
            var leaderPeptide = string.Empty;
            while (leaderboard.Any ())
            {
                var nextLeaderboard = new HashSet<string> ();
                foreach (var peptide in leaderboard)
                {
                    foreach (var component in Patterns.Peptides.Keys)
                    {
                        nextLeaderboard.Add (peptide + component);
                    }
                }
                leaderboard = nextLeaderboard;

                foreach (var peptide in leaderboard.ToArray ())
                {
                    var peptideMass = peptide.Select (p => Patterns.Peptides[p]).Sum ();
                    if (peptideMass == parentMass)
                    {
                        if (ScorePeptide (peptide, spectrum) > ScorePeptide (leaderPeptide, spectrum))
                        {
                            leaderPeptide = peptide;
                        }
                    }
                    else if (peptideMass > parentMass)
                    {
                        leaderboard.Remove (peptide);
                    }
                }

                leaderboard = new HashSet<string> (TrimLeaderboard (leaderboard, spectrum, N));
            }

            return leaderPeptide;
        }

        private static IEnumerable<string> TrimLeaderboard (HashSet<string> leaderboard, IEnumerable<int> spectrum, int N)
        {
            var scores = new Dictionary<string, int> ();
            foreach (var peptide in leaderboard)
            {
                scores[peptide] = ScorePeptide (peptide, spectrum);
            }

            var peptidesInOrder = leaderboard.OrderByDescending (p => scores[p]).ToArray ();
            var leaderboardSize = leaderboard.Count;
            if (leaderboardSize <= N)
            {
                return leaderboard;
            }

            var cutoff = scores[peptidesInOrder[N]];
            for (var j = N + 1; j < leaderboardSize; ++j)
            {
                if (scores[peptidesInOrder[j]] < cutoff)
                {
                    return peptidesInOrder.Take (j).ToArray ();
                }
            }
            return peptidesInOrder;
        }

        private static int ScorePeptide (string peptide, IEnumerable<int> spectrum)
        {
            var theoreticalSpectrum = Cyclospectrum (peptide);
            var frequency = new Dictionary<int, int> ();
            foreach (var mass in theoreticalSpectrum)
            {
                if (frequency.ContainsKey (mass))
                {
                    frequency[mass] += 1;
                }
                else
                {
                    frequency[mass] = 1;
                }
            }

            var score = 0;
            foreach (var mass in spectrum)
            {
                if (frequency.ContainsKey (mass))
                {
                    if (frequency[mass] > 0)
                    {
                        ++score;
                        --frequency[mass];
                    }
                }
            }

            return score;
        }

        private static IEnumerable<string> CyclopeptideSequencing (IEnumerable<int> spectrum)
        {
            var results = new HashSet<string> ();
            var spectrumPeptides = new HashSet<int> (spectrum);

            var parentMass = spectrum.Max ();
            var peptides = new HashSet<string> (new[] { "" });
            while (peptides.Any ())
            {
                var nextPeptides = new HashSet<string> ();
                foreach (var peptide in peptides)
                {
                    foreach (var component in Patterns.Peptides.Keys)
                    {
                        nextPeptides.Add (peptide + component);
                    }
                }
                peptides = nextPeptides;

                foreach (var peptide in peptides.ToArray ())
                {
                    var peptideMass = peptide.Select (p => Patterns.Peptides[p]).Sum ();
                    if (peptideMass == parentMass)
                    {
                        var testSpectrum = Cyclospectrum (peptide);
                        if (testSpectrum.SequenceEqual (spectrum))
                        {
                            results.Add (peptide);
                        }

                        peptides.Remove (peptide);
                    }
                    else if (!spectrumPeptides.Contains (peptideMass))
                    {
                        peptides.Remove (peptide);
                    }
                }
            }

            return results.ToArray ();
        }

        private static IEnumerable<int> Cyclospectrum (string peptide)
        {
            var prefixMass = new List<int> (new[] { 0 });

            var lastPrefixMass = 0;
            foreach (var p in peptide)
            {
                lastPrefixMass += Patterns.Peptides[p];
                prefixMass.Add (lastPrefixMass);
            }

            var peptideMass = prefixMass.Last ();
            var cyclicSpectrum = new List<int> (new[] { 0 });
            var peptideLength = peptide.Length;
            for (var i = 0; i <= peptideLength - 1; ++i)
            {
                for (var j = i + 1; j <= peptideLength; ++j)
                {
                    var mass = prefixMass[j] - prefixMass[i];
                    cyclicSpectrum.Add (mass);
                    if (i > 0 && j < peptideLength)
                    {
                        mass = peptideMass - mass;
                        cyclicSpectrum.Add (mass);
                    }
                }
            }

            return cyclicSpectrum.OrderBy (c => c).ToArray ();
        }

        private static IEnumerable<string> PeptideToRna (string peptides)
        {
            var result = new HashSet<string> (new[] { string.Empty });
            foreach (var peptide in peptides)
            {
                var nextResult = new HashSet<string> ();
                foreach (var chain in result)
                {
                    foreach (var tail in Patterns.RawCodons[peptide])
                    {
                        nextResult.Add (chain + tail);
                    }
                }
                result = nextResult;
            }
            return result;
        }

        private static IEnumerable<string> GetCodons (string dna, string peptide)
        {
            var dnaToUse = dna.Replace ('T', 'U');
            var searchQueries = new List<string> (new[] { string.Empty });
            foreach (var chunk in peptide)
            {
                var updatedQueries = new List<string> ();
                foreach (var codon in Patterns.RawCodons[chunk])
                {
                    foreach (var query in searchQueries)
                    {
                        updatedQueries.Add (query + codon);
                    }
                }
                searchQueries = updatedQueries;
            }
            searchQueries = searchQueries.Distinct ().ToList ();

            var complementQueries = new List<string> ();
            foreach (var searchQuery in searchQueries.ToArray ())
            {
                var complementQuery = ReverseComplementPattern (searchQuery.Replace ('U', 'T')).Replace ('T', 'U');
                complementQueries.Add (searchQuery);
                complementQueries.Add (complementQuery);
            }

            var result = new List<string> ();
            for (var startPosition = 0; startPosition < dna.Length; ++startPosition)
            {
                foreach (var query in complementQueries)
                {
                    if (ContainsReverseOrForward (query, dnaToUse, startPosition))
                    {
                        var pattern = dna.Substring (startPosition, query.Length);
                        result.Add (pattern);
                    }
                }
            }

            var stopSequences = Patterns.RawCodons['?'];
            return result.Where (s => !stopSequences.Any (stop => s.Contains (stop)));
        }

        private static bool ContainsReverseOrForward (string sequence, string dna, int dnaStartPosition)
        {
            var sequenceLength = sequence.Length;

            if (dnaStartPosition + sequenceLength > dna.Length)
            {
                return false;
            }

            var sequencePosition = 0;
            for (; sequencePosition < sequenceLength; ++sequencePosition)
            {
                if (sequence[sequencePosition] != dna[dnaStartPosition + sequencePosition])
                {
                    return false;
                }
            }

            return true;
        }

        private static string CircularStringProblem (int k)
        {
            var lines = new[] { k.ToString () }.Concat (GenerateWords (k)).ToArray ();
            return DeBruijnStringConstruction (lines);
        }

        private static IEnumerable<string> GenerateWords (int k)
        {
            var size = Math.Pow (2, k);
            var queue = new List<string> (new[] { "0", "1" });
            while (queue.Count != size)
            {
                var nextQueue = new List<string> (2 * queue.Count);
                foreach (var prefix in queue)
                {
                    nextQueue.Add (prefix + "0");
                    nextQueue.Add (prefix + "1");
                }
                queue = nextQueue;
            }

            return queue;
        }

        private static string PairGraphConstruction (string filepath)
        {
            var lines = File.ReadAllLines (filepath);
            var numbers = Regex.Split (lines[0], @"\s+").Select (int.Parse).ToArray ();
            var k = numbers[0];
            var d = numbers[1];

            var label = 0;
            var matrix = new Dictionary<Pair, List<Pair>> ();
            var labels = new Dictionary<int, Pair> ();
            var reverseLabels = new Dictionary<Pair, int> ();
            foreach (var line in lines.Skip (1))
            {
                var rawPair = Regex.Split (line, @"\|").ToArray ();
                var pair = new Pair (rawPair[0], rawPair[1]);

                var prefix = pair.Prefix;
                var suffix = pair.Suffix;

                if (!reverseLabels.ContainsKey (prefix))
                {
                    labels[label] = prefix;
                    reverseLabels[prefix] = label;
                    matrix[prefix] = new List<Pair> ();
                    ++label;
                }

                if (!reverseLabels.ContainsKey (suffix))
                {
                    labels[label] = suffix;
                    reverseLabels[suffix] = label;
                    matrix[suffix] = new List<Pair> ();
                    ++label;
                }

                matrix[prefix].Add (suffix);
            }

            var graph = new Graph (labels.Count);
            foreach (var row in matrix)
            {
                var source = reverseLabels[row.Key];
                foreach (var destination in row.Value.Select (dest => reverseLabels[dest]))
                {
                    graph.AddEdge (source, destination, 0);
                }
            }
            var cycle = graph.EulerianPath ();

            var length = 2 * k + d + cycle.Count () - 2;
            var path = Enumerable.Range (0, length).Select (c => new List<char> ()).ToArray ();
            var position = 0;
            foreach (var kmer in cycle)
            {
                var pair = labels[kmer];
                for (var index = 0; index < k - 1; ++index)
                {
                    var left = position + index;
                    var right = k + d + position + index;
                    path[left].Add (pair.Left[index]);
                    path[right].Add (pair.Right[index]);
                }

                ++position;
            }

            var builder = new StringBuilder ();
            foreach (var nodeProposals in path)
            {
                if (nodeProposals.Distinct ().Count () == 1)
                {
                    builder.Append (nodeProposals.First ());
                }
                else
                {
                    builder.Append ("?");
                }
            }
            return builder.ToString ();
        }

        private class Pair
        {
            public readonly string Left;
            public readonly string Right;

            public Pair (string left, string right)
            {
                Left = left;
                Right = right;
            }

            public Pair Prefix
            {
                get { return new Pair (Left.Substring (0, Left.Length - 1), Right.Substring (0, Right.Length - 1)); }
            }

            public Pair Suffix
            {
                get { return new Pair (Left.Substring (1), Right.Substring (1)); }
            }

            public bool StartsWith (Pair pair)
            {
                return Left.StartsWith (pair.Left) && Right.StartsWith (pair.Right);
            }

            public bool EndsWith (Pair pair)
            {
                return Left.EndsWith (pair.Left) && Right.EndsWith (pair.Right);
            }

            public override bool Equals (object obj)
            {
                var other = obj as Pair;
                if (other == null)
                {
                    return false;
                }
                return Left == other.Left && Right == other.Right;
            }

            public override int GetHashCode ()
            {
                return Left.GetHashCode () ^ Right.GetHashCode ();
            }
        }

        private static string DeBruijnStringConstruction (string[] lines)
        {
            var matrix = new Dictionary<string, List<string>> ();
            var labels = new Dictionary<string, int> ();
            var reverseLabels = new Dictionary<int, string> ();

            var k = int.Parse (lines.First ());
            var label = 0;
            foreach (var line in lines.Skip (1))
            {
                var prefix = line.Substring (0, k - 1);
                var suffix = line.Substring (1);

                foreach (var pattern in new[] { prefix, suffix })
                {
                    if (!labels.ContainsKey (pattern))
                    {
                        labels[pattern] = label;
                        reverseLabels[label] = pattern;
                        matrix[pattern] = new List<string> ();
                        ++label;
                    }
                }

                matrix[prefix].Add (suffix);
            }

            var graph = new Graph (label);
            foreach (var pair in matrix)
            {
                var source = labels[pair.Key];
                foreach (var destination in pair.Value.Select (key => labels[key]))
                {
                    graph.AddEdge (source, destination, 0);
                }
            }

            var path = graph.EulerianPath ().Select (lab => reverseLabels[lab]).ToArray ();
            var builder = new StringBuilder ();
            builder.Append (path.First ());
            foreach (var element in path.Skip (1))
            {
                builder.Append (element[k - 2]);
            }
            return builder.ToString ();
        }

        private static IEnumerable<string> ContingenGenerationProblem (string[] lines)
        {
            var matrix = new Dictionary<string, List<string>> ();
            var labels = new Dictionary<string, int> ();
            var reverseLabels = new Dictionary<int, string> ();

            var k = lines.First ().Length;
            var label = 0;
            foreach (var line in lines)
            {
                var prefix = line.Substring (0, k - 1);
                var suffix = line.Substring (1);

                foreach (var node in new[] { prefix, suffix })
                {
                    if (!labels.ContainsKey (node))
                    {
                        labels[node] = label;
                        reverseLabels[label] = node;
                        matrix[node] = new List<string> ();
                        ++label;
                    }
                }

                matrix[prefix].Add (suffix);
            }

            var graph = new Graph (label);
            foreach (var pair in matrix)
            {
                var source = labels[pair.Key];
                foreach (var destination in pair.Value.Select (key => labels[key]))
                {
                    graph.AddEdge (source, destination, 0);
                }
            }

            var paths = graph.MaximalNonBranchingPaths ();
            return paths.Select (path => JoinContingen (path, reverseLabels));
        }

        private static string JoinContingen (IEnumerable<int> path, Dictionary<int, string> mapping)
        {
            var builder = new StringBuilder ();
            builder.Append (mapping[path.First ()]);
            foreach (var segment in path.Skip (1))
            {
                var mappedSegment = mapping[segment];
                builder.Append (mappedSegment[mappedSegment.Length - 1]);
            }
            return builder.ToString ();
        }

        private static string EulerStringReconstruction (string[] lines)
        {
            var matrix = new Dictionary<string, List<string>> ();
            var labels = new Dictionary<string, int> ();
            var reverseLabels = new Dictionary<int, string> ();

            var k = int.Parse (lines.First ());
            var label = 0;
            foreach (var line in lines.Skip (1))
            {
                if (!labels.ContainsKey (line))
                {
                    labels[line] = label;
                    reverseLabels[label] = line;
                    matrix[line] = new List<string> ();
                }

                var prefix = line.Substring (0, k - 1);
                var suffix = line.Substring (1);
                foreach (var other in matrix.Keys)
                {
                    if (other.EndsWith (prefix))
                    {
                        matrix[other].Add (line);
                    }

                    if (other.StartsWith (suffix))
                    {
                        matrix[line].Add (other);
                    }
                }

                ++label;
            }

            var graph = new Graph (label);
            foreach (var pair in matrix)
            {
                var source = labels[pair.Key];
                foreach (var destination in pair.Value.Select (key => labels[key]))
                {
                    graph.AddEdge (source, destination, 0);
                }
            }

            var path = graph.EulerianPath ().Select (lab => reverseLabels[lab]).ToArray ();
            var builder = new StringBuilder ();
            builder.Append (path.First ());
            foreach (var element in path.Skip (1))
            {
                builder.Append (element[k - 1]);
            }
            return builder.ToString ();
        }

        private static string EulerStringReconstruction (string filepath)
        {
            var lines = File.ReadAllLines (filepath);
            return EulerStringReconstruction (lines);
        }

        private class Graph
        {
            private readonly List<Edge> edges;
            private readonly Dictionary<int, Node> nodes;

            public Graph (int nodes)
            {
                this.edges = new List<Edge> ();
                this.nodes = new Dictionary<int, Node> ();

                for (var label = 0; label < nodes; ++label)
                {
                    this.nodes[label] = new Node (label);
                }
            }

            public static Graph Parse (string path)
            {
                var lines = File.ReadAllLines (path);
                var count = 0;
                foreach (var line in lines)
                {
                    var localMax = Regex.Split (line, @"(?:\s+->\s+)|,").Select (int.Parse).Max ();
                    count = Math.Max (count, localMax + 1);
                }
                var graph = new Graph (count);

                foreach (var line in lines)
                {
                    var numbers = Regex.Split (line, @"(?:\s+->\s+)|,").Select (int.Parse).ToArray ();
                    foreach (var destination in numbers.Skip (1))
                    {
                        graph.AddEdge (numbers[0], destination, 0);
                    }
                }

                return graph;
            }

            public IEnumerable<int> EulerianPath ()
            {
                var nodeRank = new int[nodes.Count];
                foreach (var edge in edges)
                {
                    ++nodeRank[edge.Destination.Label];
                    --nodeRank[edge.Source.Label];
                }

                var endNodes = new List<Node> ();
                for (var label = 0; label < nodeRank.Length; ++label)
                {
                    if (Math.Abs (nodeRank[label]) % 2 == 1)
                    {
                        endNodes.Add (nodes[label]);
                    }
                }

                var workingCopy = Copy ();
                if (endNodes.Any ())
                {
                    var source = endNodes.Where (n => nodeRank[n.Label] > 0).First ();
                    var destination = endNodes.Where (n => nodeRank[n.Label] < 0).First ();
                    workingCopy.AddEdge (source.Label, destination.Label, 0);
                    var cycle = workingCopy.EulerianCycle ();
                    return cycle.SkipWhile (label => label != destination.Label).Concat (cycle.TakeWhile (label => label != destination.Label).Skip (1));
                }
                else
                {
                    return workingCopy.EulerianCycle ();
                }
            }

            public IEnumerable<IEnumerable<int>> MaximalNonBranchingPaths ()
            {
                var paths = new List<List<Edge>> ();

                var workingCopy = Copy ();
                var nonBranchingNodes = new HashSet<Node> (workingCopy.nodes.Values.Where (n => n.OutputEdges.Count () == 1 && n.InputEdges.Count () == 1));
                foreach (var node in workingCopy.nodes.Values)
                {
                    if (node.OutputEdges.Count () == 1 && node.InputEdges.Count () == 1)
                    {
                        continue;
                    }

                    if (node.OutputEdges.Count () > 0)
                    {
                        foreach (var sourceEdge in node.OutputEdges.ToArray ())
                        {
                            workingCopy.RemoveEdge (sourceEdge);

                            var nonBranchingPath = new List<Edge> (new[] { sourceEdge });
                            var currentEdge = sourceEdge;
                            while (nonBranchingNodes.Contains (currentEdge.Destination))
                            {
                                var edge = currentEdge.Destination.OutputEdges.First ();
                                workingCopy.RemoveEdge (edge);
                                nonBranchingPath.Add (edge);

                                currentEdge = edge;
                            }
                            paths.Add (nonBranchingPath);
                        }
                    }
                }

                while (workingCopy.edges.Any ())
                {
                    var edge = workingCopy.edges.First ();
                    var path = new List<Edge> ();
                    while (edge != null)
                    {
                        path.Add (edge);
                        workingCopy.RemoveEdge (edge);

                        edge = edge.Destination.OutputEdges.FirstOrDefault ();
                    }

                    paths.Add (path);
                }

                return paths.Select (path => new[] { path.First ().Source.Label }.Concat (path.Select (e => e.Destination.Label)));
            }

            public IEnumerable<IEnumerable<int>> Contingents ()
            {
                var workingCopy = Copy ();

                var cycles = new List<List<Edge>> ();
                var cycle = new List<Edge> ();

                Edge currentEdge = null;
                while (workingCopy.edges.Any ())
                {
                    if (cycle.Any ())
                    {
                        var node = cycle.SelectMany (e => new[] { e.Source, e.Destination }).Where (n => n.HasOutputEdges).FirstOrDefault ();
                        if (node == null)
                        {
                            node = cycles.SelectMany (c => c).SelectMany (e => new[] { e.Source, e.Destination }).Where (n => n.HasOutputEdges).FirstOrDefault ();
                        }
                        if (node == null)
                        {
                            node = workingCopy.edges.SelectMany (e => new[] { e.Source, e.Destination }).Where (n => n.HasOutputEdges).FirstOrDefault ();
                        }
                        currentEdge = node.OutputEdges.First ();
                        cycles.Add (cycle);
                    }
                    else
                    {
                        currentEdge = workingCopy.edges.FirstOrDefault ();
                    }

                    cycle = new List<Edge> ();
                    while (currentEdge != null)
                    {
                        workingCopy.RemoveEdge (currentEdge);
                        cycle.Add (currentEdge);
                        currentEdge = currentEdge.Destination.OutputEdges.FirstOrDefault ();
                    }
                }
                if (cycle.Any ())
                {
                    cycles.Add (cycle);
                }

                var pathChunks = cycles.Select (c => c.Select (edge => edge.Destination).Concat (new[] { c.First ().Destination })).Select (nodes => new PathChunk (nodes)).ToList ();
                var lastCount = pathChunks.Count + 1;
                while (pathChunks.Count < lastCount)
                {
                    lastCount = pathChunks.Count;
                    var pathToJoin = pathChunks.First ();
                    var nextChunks = new List<PathChunk> ();
                    foreach (var pathChunk in pathChunks.Skip (1))
                    {
                        if (pathToJoin.CanJoin (pathChunk))
                        {
                            pathToJoin = pathToJoin.Join (pathChunk);
                        }
                        else
                        {
                            nextChunks.Add (pathChunk);
                        }
                    }
                    nextChunks.Add (pathToJoin);
                    pathChunks = nextChunks;
                }
                return pathChunks.Select (chunk => chunk.Path).ToArray ();
            }

            public IEnumerable<int> EulerianCycle ()
            {
                var pathChunks = Contingents ();
                return pathChunks.First ();
            }

            private class PathChunk
            {
                private readonly List<Node> chunk;
                private readonly HashSet<Node> nodes;

                public PathChunk (IEnumerable<Node> nodes)
                {
                    this.chunk = new List<Node> (nodes);
                    this.nodes = new HashSet<Node> (nodes);
                }

                public bool CanJoin (PathChunk other)
                {
                    return nodes.Intersect (other.nodes).Any ();
                }

                public IEnumerable<int> Path
                {
                    get { return chunk.Select (n => n.Label); }
                }

                public PathChunk Join (PathChunk other)
                {
                    var commonNode = nodes.Intersect (other.nodes).First ();
                    var firstSegment = chunk.TakeWhile (node => node != commonNode).ToArray ();
                    var secondSegment = other.chunk.SkipWhile (node => node != commonNode).ToArray ();
                    var thirdSegment = other.chunk.Skip (1).TakeWhile (node => node != commonNode);
                    var lastSegment = chunk.SkipWhile (node => node != commonNode);
                    return new PathChunk (firstSegment.Concat (secondSegment).Concat (thirdSegment).Concat (lastSegment));
                }
            }

            public void AddEdge (int source, int destination, int weight)
            {
                var sourceNode = nodes[source];
                var destinationNode = nodes[destination];
                var edge = new Edge (sourceNode, destinationNode, weight);
                edges.Add (edge);
                sourceNode.AddOutputEdge (edge);
                destinationNode.AddInputEdge (edge);
            }

            public void RemoveEdge (Edge edge)
            {
                this.edges.Remove (edge);

                var sourceNode = nodes[edge.Source.Label];
                var destinationNode = nodes[edge.Destination.Label];
                sourceNode.RemoveEdge (edge);
                destinationNode.RemoveEdge (edge);
            }

            public Graph Copy ()
            {
                var copy = new Graph (nodes.Count);
                foreach (var edge in edges)
                {
                    copy.AddEdge (edge.Source.Label, edge.Destination.Label, edge.Weight);
                }
                return copy;
            }

            public Tuple<IEnumerable<int>, int> LongestPath ()
            {
                var parentDict = new Dictionary<Node, Node> ();
                var weightDict = new Dictionary<Node, int> ();
                foreach (var node in nodes.Values)
                {
                    weightDict[node] = 0;
                }

                foreach (var node in Sort ())
                {
                    foreach (var edge in node.OutputEdges)
                    {
                        var candidateWeight = weightDict[node] + edge.Weight;
                        if (weightDict[edge.Destination] <= candidateWeight)
                        {
                            weightDict[edge.Destination] = candidateWeight;
                            parentDict[edge.Destination] = node;
                        }
                    }
                }

                var max = weightDict.Values.Max ();
                var lastNode = weightDict.Where (pair => pair.Value == max).First ().Key;
                var path = new List<Node> ();
                var currentNode = lastNode;
                while (currentNode != null)
                {
                    path.Add (currentNode);
                    parentDict.TryGetValue (currentNode, out currentNode);
                }

                return new Tuple<IEnumerable<int>, int> (path.Select (n => n.Label).Reverse ().ToArray (), weightDict[lastNode]);
            }

            private IEnumerable<Node> Sort ()
            {
                var workingCopy = Copy ();
                var result = new List<Node> ();
                var sourceNodes = new HashSet<Node> ();
                foreach (var node in workingCopy.nodes.Values.Where (n => !n.HasInputEdges))
                {
                    sourceNodes.Add (node);
                }

                while (sourceNodes.Any ())
                {
                    var sourceNode = sourceNodes.First ();
                    sourceNodes.Remove (sourceNode);

                    result.Add (sourceNode);
                    foreach (var edge in sourceNode.OutputEdges.ToArray ())
                    {
                        workingCopy.RemoveEdge (edge);

                        if (!edge.Destination.HasInputEdges)
                        {
                            sourceNodes.Add (edge.Destination);
                        }
                    }
                }

                if (workingCopy.edges.Any ())
                {
                    throw new InvalidOperationException ();
                }

                return result.Select (n => nodes[n.Label]).ToArray ();
            }
        }

        private class Node
        {
            private readonly int label;

            private readonly List<Edge> inputEdges;
            private readonly List<Edge> outputEdges;

            public Node (int label)
            {
                this.label = label;

                this.inputEdges = new List<Edge> ();
                this.outputEdges = new List<Edge> ();
            }

            public IEnumerable<Edge> InputEdges { get { return inputEdges; } }

            public IEnumerable<Edge> OutputEdges { get { return outputEdges; } }

            public int Label { get { return label; } }

            public bool HasInputEdges { get { return inputEdges.Any (); } }

            public bool HasOutputEdges { get { return outputEdges.Any (); } }

            public bool RemoveEdge (Edge edge)
            {
                var result = false;
                result |= inputEdges.Remove (edge);
                result |= outputEdges.Remove (edge);
                return result;
            }

            public void AddInputEdge (Edge edge)
            {
                inputEdges.Add (edge);
            }

            public void AddOutputEdge (Edge edge)
            {
                outputEdges.Add (edge);
            }

            public override bool Equals (object obj)
            {
                var other = obj as Node;
                return label == other.label;
            }

            public override int GetHashCode ()
            {
                return label;
            }

            public override string ToString ()
            {
                return label.ToString ();
            }
        }

        private class Edge
        {
            private readonly Node source;
            private readonly Node destination;
            private readonly int weight;

            public Edge (Node source, Node destination, int weight)
            {
                this.source = source;
                this.destination = destination;
                this.weight = weight;
            }

            public Node Source { get { return source; } }
            public Node Destination { get { return destination; } }
            public int Weight { get { return weight; } }

            public override bool Equals (object obj)
            {
                var other = obj as Edge;
                if (other == null)
                {
                    return false;
                }

                return Source == other.Source
                    && Destination == other.Destination
                    && Weight == other.Weight;
            }

            public override int GetHashCode ()
            {
                return source.GetHashCode () ^ destination.GetHashCode () ^ weight;
            }
        }

        private class ManhattanTourist
        {
            public static int Solve (string filepath)
            {
                var lines = File.ReadAllLines (filepath);
                var rawCoordinates = lines.First ();
                var coordinaes = Regex.Split (rawCoordinates, @"\s+").Select (int.Parse).ToArray ();
                var rows = coordinaes[0];
                var columns = coordinaes[1];

                var downMatrix = lines.Skip (1).Take (rows).Select (r => Regex.Split (r, @"\s+").Select (int.Parse).ToArray ()).ToArray ();
                var rightMatrix = lines.Skip (2 + rows).Select (r => Regex.Split (r, @"\s+").Select (int.Parse).ToArray ()).ToArray ();

                var solutionMatrix = new int[rows + 1, columns + 1];
                solutionMatrix[0, 0] = 0;

                for (var i = 1; i <= rows; ++i)
                {
                    solutionMatrix[i, 0] = solutionMatrix[i - 1, 0] + downMatrix[i - 1][0];
                }

                for (var j = 1; j <= columns; ++j)
                {
                    solutionMatrix[0, j] = solutionMatrix[0, j - 1] + rightMatrix[0][j - 1];
                }

                for (var i = 1; i <= rows; ++i)
                {
                    for (var j = 1; j <= columns; ++j)
                    {
                        solutionMatrix[i, j] = Math.Max (solutionMatrix[i - 1, j] + downMatrix[i - 1][j], solutionMatrix[i, j - 1] + rightMatrix[i][j - 1]);
                    }
                }

                return solutionMatrix[rows, columns];
            }
        }

        private class LCSBackTrack
        {
            private enum Direction { South, East, SouthEast }

            public static string Compare (string left, string right)
            {
                var backtrack = BackTrack (left, right);
                return Output (backtrack, left);
            }

            private static string Output (Direction[,] backtrack, string left)
            {
                var builder = new StringBuilder ();

                var i = backtrack.GetLength (0) - 1;
                var j = backtrack.GetLength (1) - 1;
                while (i != 0 && j != 0)
                {
                    switch (backtrack[i, j])
                    {
                        case Direction.SouthEast:
                            builder.Append (left[i - 1]);
                            --i;
                            --j;
                            break;
                        case Direction.East:
                            --j;
                            break;
                        case Direction.South:
                            --i;
                            break;
                    }
                }

                return string.Join (string.Empty, builder.ToString ().Reverse ());
            }

            private static Direction[,] BackTrack (string left, string right)
            {
                var v = left.Length + 1;
                var w = right.Length + 1;
                var s = new int[v, w];
                var backtrack = new Direction[v, w];
                for (var i = 0; i < v; ++i)
                {
                    s[i, 0] = 0;
                }
                for (var j = 0; j < w; ++j)
                {
                    s[0, j] = 0;
                }

                for (var i = 1; i < v; ++i)
                {
                    for (var j = 1; j < w; ++j)
                    {
                        s[i, j] = Math.Max (s[i - 1, j], s[i, j - 1]);
                        if (left[i - 1] == right[j - 1])
                        {
                            s[i, j] = Math.Max (s[i, j], s[i - 1, j - 1] + 1);
                        }

                        if (s[i, j] == s[i - 1, j])
                        {
                            backtrack[i, j] = Direction.South;
                        }
                        else if (s[i, j] == s[i, j - 1])
                        {
                            backtrack[i, j] = Direction.East;
                        }
                        else if (s[i, j] == (s[i - 1, j - 1] + 1) && left[i - 1] == right[j - 1])
                        {
                            backtrack[i, j] = Direction.SouthEast;
                        }
                    }
                }

                return backtrack;
            }
        }

        private class DPChange
        {
            private readonly int[] coins;
            private readonly Dictionary<int, int> change;

            public DPChange (IEnumerable<int> coins)
            {
                this.coins = coins.ToArray ();
                this.change = new Dictionary<int, int> ();
            }

            public int Change (int money)
            {
                if (money == 0)
                {
                    return 0;
                }

                var count = 0;
                if (change.TryGetValue (money, out count))
                {
                    return count;
                }

                count = int.MaxValue;
                foreach (var coin in coins)
                {
                    if (money >= coin)
                    {
                        var localCount = Change (money - coin) + 1;
                        if (localCount < count)
                        {
                            count = localCount;
                        }
                    }
                }

                var savedCount = int.MaxValue;
                if (change.TryGetValue (money, out savedCount))
                {
                    if (count < savedCount)
                    {
                        change[money] = count;
                    }
                }
                else if (count != int.MaxValue)
                {
                    change[money] = count;
                }

                return count;
            }
        }

        private static IEnumerable<string[]> AdjacencyLists (string[] dna)
        {
            var chunks = dna.Select (d => new Chunk (d)).ToArray ();
            var result = new List<string[]> ();
            for (var position = 0; position < chunks.Length; ++position)
            {
                var previous = chunks[position];
                var adjacencyList = new List<string> ();
                adjacencyList.Add (previous.Data);
                foreach (var next in Enumerable.Concat (chunks.Take (position - 1), chunks.Skip (position)))
                {
                    if (next.Follows (previous))
                    {
                        adjacencyList.Add (next.Data);
                    }
                }

                if (adjacencyList.Count > 1)
                {
                    result.Add (adjacencyList.ToArray ());
                }
            }
            return result;
        }

        private static IEnumerable<string> DeBruijnFormat (string[] kmers)
        {
            var k = kmers.First ().Length;
            var matrix = new Dictionary<string, List<string>> ();
            foreach (var kmer in kmers)
            {
                var prev = kmer.Substring (0, k - 1);
                var next = kmer.Substring (1);

                List<string> rowToUse = null;
                if (!matrix.ContainsKey (prev))
                {
                    rowToUse = matrix[prev] = new List<string> ();
                }
                else
                {
                    rowToUse = matrix[prev];
                }

                rowToUse.Add (next);
            }

            var result = new List<string> ();
            foreach (var pair in matrix)
            {
                var neighbours = pair.Value.OrderBy (i => i).ToArray ();
                result.Add (pair.Key + " -> " + string.Join (",", neighbours));
            }
            return result.OrderBy (r => r);
        }

        private static IEnumerable<string> DeBruijnFormat (string dna, int k)
        {
            var kmers = EnumerateKMers (dna, k).ToArray ();
            return DeBruijnFormat (kmers);
        }

        private static string NaiveConstruction (string[] chunks)
        {
            var builder = new StringBuilder ();
            builder.Append (chunks[0]);
            for (var position = 1; position < chunks.Length; ++position)
            {
                var chunk = chunks[position];
                builder.Append (chunk.Last ());
            }
            return builder.ToString ();
        }

        private class StringReconstruction
        {
            public static string Reconstruct (Chunk[] chunks)
            {
                var acc = new List<Chunk> (chunks);
                var candidates = new HashSet<Chunk> (chunks);

                foreach (var candidate in candidates)
                {
                    acc.Remove (candidate);

                    var result = ReconstructInternal (candidate, acc);
                    if (result.Any ())
                    {
                        var builder = new StringBuilder ();
                        builder.Append (candidate.Data);
                        foreach (var chunk in result)
                        {
                            builder.Append (chunk.Last);
                        }
                        return builder.ToString ();
                    }

                    acc.Add (candidate);
                }

                return string.Empty;
            }

            public static LinkedList<Chunk> ReconstructInternal (Chunk candidate, List<Chunk> acc)
            {
                if (!acc.Any ())
                {
                    var result = new LinkedList<Chunk> ();
                    result.AddLast (candidate);
                    return result;
                }

                var follows = acc.Where (c => c.Follows (candidate)).Distinct ().ToArray ();
                foreach (var follow in follows)
                {
                    acc.Remove (follow);
                    var result = ReconstructInternal (follow, acc);
                    if (result.Any ())
                    {
                        result.AddFirst (candidate);
                        return result;
                    }
                    acc.Add (follow);
                }

                return new LinkedList<Chunk> ();
            }
        }

        private class Chunk
        {
            private readonly string _chunk;

            public Chunk (string chunk)
            {
                _chunk = chunk;
            }

            public override bool Equals (object obj)
            {
                var chunk = obj as Chunk;
                if (chunk == null)
                {
                    return false;
                }

                if (Length != chunk.Length)
                {
                    return false;
                }

                var size = Length;
                for (var position = 0; position < size; ++position)
                {
                    if (_chunk[position] != chunk._chunk[position])
                    {
                        return false;
                    }
                }
                return true;
            }

            public override int GetHashCode ()
            {
                var result = 0;
                foreach (var nucleotide in _chunk)
                {
                    result ^= nucleotide;
                }
                return result;
            }

            public bool Follows (Chunk previous)
            {
                var size = Math.Min (Length, previous.Length);
                for (int previousPos = 1, nextPos = 0; previousPos < size; ++previousPos, ++nextPos)
                {
                    if (_chunk[nextPos] != previous._chunk[previousPos])
                    {
                        return false;
                    }
                }
                return true;
            }

            public char this[int index]
            {
                get
                {
                    return _chunk[index];
                }
            }

            public char Last { get { return _chunk[Length - 1]; } }

            public int Length { get { return _chunk.Length; } }

            public string Data { get { return _chunk; } }
        }

        private static string[] GibbsSampler (string[] dna, int k, int t, int N, int times)
        {
            var score = int.MaxValue;
            var bestMotif = new string[0];

            for (var iteration = 0; iteration < times; ++iteration)
            {
                var localMotif = GibbsSampler (dna, k, t, N);
                var localScore = Score (localMotif);

                if (localScore < score)
                {
                    score = localScore;
                    bestMotif = localMotif;
                }
            }

            return bestMotif;
        }

        private static string[] GibbsSampler (string[] dna, int k, int t, int N)
        {
            var size = dna[0].Length - k;
            var random = new Random ((int)DateTime.Now.Ticks);
            var motifs = new List<string> ();
            foreach (var row in dna)
            {
                var start = random.Next (size);
                motifs.Add (row.Substring (start, k));
            }
            var bestMotifs = motifs.ToArray ();

            for (var j = 1; j < N; ++j)
            {
                var i = random.Next (t);

                var position = 0;
                var profileBuilder = new ProfileBuilder ();
                foreach (var motif in bestMotifs)
                {
                    if (position != i)
                    {
                        profileBuilder.Add (motif);
                    }

                    ++position;
                }
                var profile = profileBuilder.Build ();

                motifs = bestMotifs.Take (i - 1).ToList ();
                motifs.Add (profile.GenerateKMer (dna[i]));
                motifs.AddRange (bestMotifs.Skip (i));

                var localMotifs = motifs.ToArray ();
                if (Score (localMotifs) < Score (bestMotifs))
                {
                    bestMotifs = localMotifs;
                }
            }

            return bestMotifs;
        }

        private static string[] RandomizedMotifSearch (string[] dna, int k, int times)
        {
            var score = int.MaxValue;
            var bestMotif = new string[0];

            for (var iteration = 0; iteration < times; ++iteration)
            {
                var localMotif = RandomizedMotifSearch (dna, k);
                var localScore = Score (localMotif);

                if (localScore < score)
                {
                    score = localScore;
                    bestMotif = localMotif;
                }
            }

            return bestMotif;
        }

        private static string[] RandomizedMotifSearch (string[] dna, int k)
        {
            var size = dna[0].Length - k;
            var random = new Random ((int)DateTime.Now.Ticks);
            var motifs = new List<string> ();
            /*foreach (var row in dna)
            {
                var start = random.Next (size);
                motifs.Add (row.Substring (start, k));
            }
            var bestMotifs = motifs.ToArray ();*/
            var bestMotifs = new[] {
                "GTC",
                "CCC",
                "ATA",
                "GCT"
            };
            motifs = bestMotifs.ToList ();

            while (true)
            {
                var profileBuilder = new ProfileBuilder ();
                motifs.ForEach (m => profileBuilder.Add (m));
                var profile = profileBuilder.Build ();

                motifs = new List<string> ();
                foreach (var row in dna)
                {
                    motifs.Add (profile.MostProbableKmer (row));
                }

                var localMotifs = motifs.ToArray ();
                if (Score (localMotifs) < Score (bestMotifs))
                {
                    bestMotifs = localMotifs;
                }
                else
                {
                    return bestMotifs;
                }
            }
        }

        private static string[] GreedyMotifSearch (string[] dna, int k, int t)
        {
            var bestMotifs = dna.Select (line => line.Substring (0, k)).ToArray ();
            foreach (var kmer in EnumerateKMers (dna[0], k))
            {
                var candidateMotifs = new string[dna.Length];
                candidateMotifs[0] = kmer;
                for (var i = 1; i < t; ++i)
                {
                    var profileBuilder = new ProfileBuilder ();
                    for (var motifIndex = 0; motifIndex < i; ++motifIndex)
                    {
                        profileBuilder.Add (candidateMotifs[motifIndex]);
                    }
                    var profile = profileBuilder.Build ();
                    candidateMotifs[i] = profile.MostProbableKmer (dna[i]);
                }

                if (Score (candidateMotifs) < Score (bestMotifs))
                {
                    bestMotifs = candidateMotifs;
                }
            }
            return bestMotifs;
        }

        private static int Score (string[] motifs)
        {
            var score = 0;
            var columns = motifs[0].Length;
            var counts = new int[4];

            for (var column = 0; column < columns; ++column)
            {
                for (var i = 0; i < 4; ++i)
                {
                    counts[i] = 0;
                }

                foreach (var row in motifs)
                {
                    var index = Patterns.ToNumber (row[column]);
                    ++counts[index];
                }

                score += counts.Sum () - counts.Max ();
            }

            return score;
        }

        private class ProfileBuilder
        {
            private List<string> _motifs;

            public ProfileBuilder ()
            {
                _motifs = new List<string> ();
            }

            public ProfileBuilder Add (string motif)
            {
                _motifs.Add (motif);

                return this;
            }

            public Profile Build ()
            {
                var columns = new List<double[]> ();
                var totalMotifs = _motifs.Count;
                var rowLength = _motifs.First ().Length;

                for (var position = 0; position < rowLength; ++position)
                {
                    var counts = new double[4];
                    foreach (var motif in _motifs)
                    {
                        var index = Patterns.ToNumber (motif[position]);
                        counts[index] += 1.0;
                    }
                    for (var i = 0; i < 4; ++i)
                    {
                        counts[i] += 1;
                    }

                    columns.Add (counts.Select (c => c / totalMotifs).ToArray ());
                }

                var rows = new List<double[]> ();
                for (var row = 0; row < 4; ++row)
                {
                    var rowContent = new double[rowLength];
                    for (var column = 0; column < rowLength; ++column)
                    {
                        rowContent[column] = columns[column][row];
                    }
                    rows.Add (rowContent);
                }

                return new Profile (rows.ToArray ());
            }
        }

        private class Profile
        {
            private readonly double[][] _matrix;

            public Profile (double[][] matrix)
            {
                _matrix = matrix;
            }

            public string Signature
            {
                get
                {
                    var index = -1;

                    var rowSize = RowSize;
                    var columnSize = ColumnSize;
                    var kmer = new StringBuilder ();
                    for (var column = 0; column < columnSize; ++column)
                    {
                        for (var row = 0; row < rowSize; ++row)
                        {
                            var probability = 0.0;
                            if (_matrix[row][column] > probability)
                            {
                                probability = _matrix[row][column];
                                index = row;
                            }
                        }
                        kmer.Append (Patterns.ToSymbol (index));
                    }

                    return kmer.ToString ();
                }
            }

            public double Pr (string text)
            {
                var result = 1.0;
                for (var column = 0; column < text.Length; ++column)
                {
                    var number = Patterns.ToNumber (text[column]);
                    result *= _matrix[number][column];
                }
                return result;
            }

            public bool IsConsensus (string text)
            {
                var characters = text.ToCharArray ();
                for (var column = 0; column < characters.Length; ++column)
                {
                    var max = MaxProbability (column);
                    var index = Patterns.ToNumber (characters[column]);
                    if (_matrix[index][column] != max)
                    {
                        return false;
                    }
                }
                return true;
            }

            private double MaxProbability (int column)
            {
                var result = 0.0;
                foreach (var row in _matrix)
                {
                    result = Math.Max (result, row[column]);
                }
                return result;
            }

            public string MostProbableKmer (string text)
            {
                var probability = 0.0;
                string mostProbableKmer = null;
                foreach (var kmer in EnumerateKMers (text, K))
                {
                    var candidateProbability = GetProbability (kmer);
                    if (candidateProbability > probability)
                    {
                        probability = candidateProbability;
                        mostProbableKmer = kmer;
                    }
                }

                if (mostProbableKmer != null)
                {
                    return mostProbableKmer;
                }

                return text.Substring (0, K);
            }

            public string GenerateKMer (string dna)
            {
                var scores = new List<double> ();
                var size = dna.Length - K + 1;
                for (var start = 0; start < size; ++start)
                {
                    var kmer = dna.Substring (start, K);
                    scores.Add (Pr (kmer));
                }

                var random = new Random ((int)DateTime.UtcNow.Ticks);
                var choice = random.NextDouble () * scores.Sum ();
                var position = 0;
                var threshold = 0.0;
                foreach (var score in scores)
                {
                    threshold += score;

                    if (choice <= threshold)
                    {
                        break;
                    }

                    ++position;
                }

                return dna.Substring (position, K);
            }

            private double GetProbability (string kmer)
            {
                var probability = 1.0;
                var position = 0;
                foreach (var nucleotide in kmer)
                {
                    var index = Patterns.ToNumber (nucleotide);
                    probability *= _matrix[index][position];
                    ++position;
                }
                return probability;
            }

            private int RowSize { get { return _matrix.Length; } }

            private int ColumnSize { get { return _matrix[0].Length; } }

            public int K { get { return ColumnSize; } }
        }

        private static IEnumerable<string> MedianStrings (string[] dna, int k)
        {
            var distance = int.MaxValue;
            var size = Math.Pow (4, k);
            var medians = new HashSet<string> ();
            for (var number = 0; number < size; ++number)
            {
                var pattern = Patterns.NumberToPattern (number, k);
                var localDistance = DistanceBetweenPatternAndStrings (pattern, dna);
                if (localDistance < distance)
                {
                    distance = localDistance;
                    medians.Clear ();
                }

                if (localDistance == distance)
                {
                    medians.Add (pattern);
                }
            }
            return medians;
        }

        private static string MedianString (string[] dna, int k)
        {
            var distance = int.MaxValue;
            var size = Math.Pow (4, k);
            var median = string.Empty;
            for (var number = 0; number < size; ++number)
            {
                var pattern = Patterns.NumberToPattern (number, k);
                var localDistance = DistanceBetweenPatternAndStrings (pattern, dna);
                if (localDistance < distance)
                {
                    distance = localDistance;
                    median = pattern;
                }
            }
            return median;
        }

        private static int DistanceBetweenPatternAndStrings (string pattern, string[] dna)
        {
            var k = pattern.Length;
            var distance = 0;
            foreach (var text in dna)
            {
                var minHammingDistance = int.MaxValue;
                foreach (var patternCandidate in EnumerateKMers (text, k))
                {
                    minHammingDistance = Math.Min (minHammingDistance, HammingDistance (patternCandidate, pattern));
                }
                distance += minHammingDistance;
            }
            return distance;
        }

        private static IEnumerable<string> EnumerateKMers (string text, int k)
        {
            var size = text.Length - k;
            for (var start = 0; start <= size; ++start)
            {
                yield return text.Substring (start, k);
            }
        }

        private static IEnumerable<string> MotifEnumeration (string[] dnaStrings, int k, int d)
        {
            var patterns = new HashSet<string> ();
            var dna = string.Join (string.Empty, dnaStrings);

            foreach (var kmer in GetKMers (dna, k))
            {
                foreach (var neighbour in Neighbours (kmer, d))
                {
                    var foundInEveryDnaString = dnaStrings.All (s => ComputeCount (neighbour, s, d) > 0);
                    if (foundInEveryDnaString)
                    {
                        patterns.Add (neighbour);
                    }
                }
            }

            return patterns;
        }

        private static IEnumerable<string> GetKMers (string dna, int k)
        {
            var size = dna.Length - k;
            for (var start = 0; start < size; ++start)
            {
                yield return dna.Substring (start, k);
            }
        }

        private static IEnumerable<string> FrequentWordsWithMismatches (string text, int k, int d)
        {
            var neighborHoods = new List<string> ();
            for (var i = 0; i < text.Length - k; ++i)
            {
                foreach (var neighbor in IterativeNeighbors (text.Substring (i, k), d))
                {
                    neighborHoods.Add (neighbor);
                }
            }

            var neighborHoodArray = neighborHoods.ToArray ();
            var count = new int[neighborHoodArray.Length];
            var sortedIndex = new int[neighborHoodArray.Length];
            for (var i = 0; i < neighborHoodArray.Length; ++i)
            {
                var pattern = neighborHoodArray[i];
                sortedIndex[i] = Patterns.ToNumber (pattern);
                count[i] = 1;
            }
            Array.Sort (sortedIndex);

            for (var i = 0; i < neighborHoodArray.Length - 1; ++i)
            {
                if (sortedIndex[i] == sortedIndex[i + 1])
                {
                    count[i + 1] = count[i] + 1;
                }
            }

            var frequentPatterns = new HashSet<string> ();
            var maxCount = count.Max ();
            for (var i = 0; i < neighborHoodArray.Length - 1; ++i)
            {
                if (count[i] == maxCount)
                {
                    var pattern = Patterns.NumberToPattern (sortedIndex[i], k);
                    frequentPatterns.Add (pattern);
                }
            }
            return frequentPatterns.ToArray ();
        }

        private static IEnumerable<string> FrequentWordsWithMismatchesAndReverseComponents (string text, int k, int d)
        {
            var neighborHoods = new List<string> ();
            for (var i = 0; i < text.Length - k; ++i)
            {
                foreach (var neighbor in IterativeNeighbors (text.Substring (i, k), d))
                {
                    neighborHoods.Add (neighbor);
                }
            }

            var neighborHoodArray = neighborHoods.ToArray ();
            var count = new int[neighborHoodArray.Length];
            var sortedIndex = new int[neighborHoodArray.Length];
            for (var i = 0; i < neighborHoodArray.Length; ++i)
            {
                var pattern = neighborHoodArray[i];
                sortedIndex[i] = Patterns.ToNumber (pattern);
                count[i] = 1;
            }
            Array.Sort (sortedIndex);

            var numberOccurence = new Dictionary<int, int> ();
            var lastIndex = neighborHoodArray.Length - 1;
            for (var i = 0; i < lastIndex; ++i)
            {
                if (sortedIndex[i] == sortedIndex[i + 1])
                {
                    count[i + 1] = count[i] + 1;
                }
                else
                {
                    numberOccurence[sortedIndex[i]] = count[i];
                }
            }
            numberOccurence[sortedIndex[lastIndex]] = count[lastIndex];

            var totalCount = new int[count.Length];
            var frequentPatterns = new HashSet<string> ();
            for (var i = 0; i < lastIndex; ++i)
            {
                if (sortedIndex[i] >= sortedIndex[i + 1])
                {
                    var pattern = Patterns.NumberToPattern (sortedIndex[i], k);
                    var rcPattern = ReverseComplementPattern (pattern);
                    var rcIndex = Patterns.ToNumber (rcPattern);
                    totalCount[i] = count[i];
                    if (numberOccurence.ContainsKey (rcIndex))
                    {
                        totalCount[i] += numberOccurence[rcIndex];
                    }
                }
            }
            var lastPattern = Patterns.NumberToPattern (sortedIndex[lastIndex], k);
            var lastRcPattern = ReverseComplementPattern (lastPattern);
            var lastRcIndex = Patterns.ToNumber (lastRcPattern);
            totalCount[lastIndex] = count[lastIndex];
            if (numberOccurence.ContainsKey (lastRcIndex))
            {
                totalCount[lastIndex] += numberOccurence[lastRcIndex];
            }

            var max = totalCount.Max ();
            for (var i = 0; i < neighborHoodArray.Length; ++i)
            {
                if (totalCount[i] == max)
                {
                    var pattern = Patterns.NumberToPattern (sortedIndex[i], k);
                    frequentPatterns.Add (pattern);

                    var rcPattern = ReverseComplementPattern (pattern);
                    var rcIndex = Patterns.ToNumber (rcPattern);
                    if (numberOccurence.ContainsKey (rcIndex))
                    {
                        frequentPatterns.Add (rcPattern);
                    }
                }
            }

            return frequentPatterns.ToArray ();
        }

        private static double Enthropy (string[] dna)
        {
            var columns = dna[0].Length;
            var rows = dna.Length;
            var columnEntropies = new double[columns];
            var counter = new int[4];

            for (var column = 0; column < columns; ++column)
            {
                for (var position = 0; position < counter.Length; ++position)
                {
                    counter[position] = 0;
                }

                foreach (var dnaString in dna)
                {
                    ++counter[Patterns.ToNumber (dnaString[column])];
                }

                columnEntropies[column] = -1.0 * counter.Select (c => GetEnthropyComponent ((double)c / rows)).Sum ();
            }

            return columnEntropies.Sum ();
        }

        private static double GetEnthropyComponent (double p)
        {
            if (p == 0.0)
            {
                return 0.0;
            }
            else
            {
                return p * Math.Log (p, 2.0);
            }
        }

        private static int ComputeCount (string pattern, string text, int d)
        {
            return ApproximatePatternMatchingProblem (pattern, text, d).Count ();
        }

        private static IEnumerable<int> ApproximatePatternMatchingProblem (string pattern, string text, int d)
        {
            var result = new List<int> ();
            var rawText = text.ToCharArray ();
            var size = rawText.Length - pattern.Length;
            for (var start = 0; start <= size; ++start)
            {
                if (HammingDistance (pattern, rawText, start) <= d)
                {
                    result.Add (start);
                }
            }
            return result;
        }

        private static int HammingDistance (string pattern, char[] text, int start)
        {
            var distance = 0;
            for (var position = 0; position < pattern.Length; ++position)
            {
                if (pattern[position] != text[start + position])
                {
                    ++distance;
                }
            }
            return distance;
        }

        private static int HammingDistance (IEnumerable<char> pattern, string text)
        {
            var distance = 0;
            var position = 0;
            foreach (var chunk in pattern)
            {
                if (chunk != text[position])
                {
                    ++distance;
                }
                ++position;
            }
            return distance;
        }

        private static int HammingDistance (string left, string right)
        {
            if (left.Length != right.Length)
            {
                return int.MaxValue;
            }

            var distance = 0;
            for (var position = 0; position < left.Length; ++position)
            {
                if (left[position] != right[position])
                {
                    ++distance;
                }
            }
            return distance;
        }

        private static IEnumerable<int> MinimumSkew (string chain)
        {
            var skew = Skew (chain);

            var result = new List<int> ();
            var min = skew.Min ();
            var position = 0;
            foreach (var occurence in skew)
            {
                if (occurence == min)
                {
                    result.Add (position);
                }
                ++position;
            }
            return result;
        }

        private static IEnumerable<int> Skew (string chain)
        {
            var counter = 0;
            var result = new List<int> ();

            foreach (var chunk in chain)
            {
                result.Add (counter);

                if (chunk == 'C')
                {
                    --counter;
                }
                else if (chunk == 'G')
                {
                    ++counter;
                }
            }
            result.Add (counter);

            return result;
        }

        private static IEnumerable<string> IterativeNeighbors (string pattern, int d)
        {
            var result = new HashSet<string> ();
            result.Add (pattern);

            for (var j = 0; j < d; ++j)
            {
                var immediateNeighbors = new HashSet<string> ();
                foreach (var neighbor in result)
                {
                    immediateNeighbors.UnionWith (ImmediateNeighbors (neighbor));
                }
                result.UnionWith (immediateNeighbors);
            }

            return result.ToArray ();
        }

        private static IEnumerable<string> ImmediateNeighbors (string pattern)
        {
            if (pattern.Count () == 1)
            {
                return Patterns.Nucleotides;
            }

            var neighbors = new HashSet<string> ();
            neighbors.Add (pattern);

            var size = pattern.Count () - 1;
            foreach (var nucleotide in Patterns.FilterNucleotides (pattern[0]))
            {
                var neighbor = nucleotide + pattern.Substring (1);
                neighbors.Add (neighbor);
            }

            for (var position = 1; position < size; ++position)
            {
                var symbol = pattern[position];
                foreach (var nucleotide in Patterns.FilterNucleotides (symbol))
                {
                    var neighbor = pattern.Substring (0, position) + nucleotide + pattern.Substring (position + 1);
                    neighbors.Add (neighbor);
                }
            }

            foreach (var nucleotide in Patterns.FilterNucleotides (pattern[size]))
            {
                var neighbor = pattern.Substring (0, size) + nucleotide;
                neighbors.Add (neighbor);
            }

            return neighbors.ToArray ();
        }

        private static IEnumerable<string> Neighbours (IEnumerable<char> pattern, int d)
        {
            if (d == 0)
            {
                return new[] { string.Join (string.Empty, pattern) };
            }

            if (pattern.Count () == 1)
            {
                return Patterns.Nucleotides;
            }

            var neighbors = new HashSet<string> ();
            var suffixNeighbors = Neighbours (pattern.Skip (1), d);
            foreach (var text in suffixNeighbors)
            {
                if (HammingDistance (pattern.Skip (1), text) < d)
                {
                    foreach (var nucleotide in Patterns.Nucleotides)
                    {
                        neighbors.Add (nucleotide + text);
                    }
                }
                else
                {
                    neighbors.Add (pattern.First () + text);
                }
            }
            return neighbors.ToArray ();
        }

        private static IEnumerable<string> ClumpFindingProblem (string genome, int k, int t, int L)
        {
            var spaceSize = (int)Math.Pow (4, k);
            var clump = new bool[spaceSize];
            var text = genome.Substring (0, L);
            var frequencyArray = ComputeFrequencies (text, k);

            for (int i = 0; i < spaceSize; ++i)
            {
                if (frequencyArray[i] >= t)
                {
                    clump[i] = true;
                }
            }

            var size = genome.Length - L + 1;
            for (int i = 1; i < size; ++i)
            {
                var firstPattern = genome.Substring (i - 1, k);
                var index = Patterns.ToNumber (firstPattern);
                frequencyArray[index] -= 1;

                var lastPattern = genome.Substring (i + L - k, k);
                index = Patterns.ToNumber (lastPattern);
                frequencyArray[index] += 1;

                if (frequencyArray[index] >= t)
                {
                    clump[index] = true;
                }
            }

            var frequentPatterns = new HashSet<string> ();
            for (int i = 0; i < spaceSize; ++i)
            {
                if (clump[i])
                {
                    var pattern = Patterns.NumberToPattern (i, k);
                    frequentPatterns.Add (pattern);
                }
            }
            return frequentPatterns.ToArray ().OrderBy (i => i).ToArray ();
        }

        private static int[] ComputeFrequencies (string text, int k)
        {
            var frequencyArray = new int[(int)Math.Pow (4, k)];
            var size = text.Length - k + 1;
            for (var i = 0; i < size; ++i)
            {
                var pattern = text.Substring (i, k);
                var j = Patterns.ToNumber (pattern);
                frequencyArray[j] += 1;
            }
            return frequencyArray;
        }

        private static string ReverseComplementPattern (string pattern)
        {
            var rcPattern = pattern.Reverse ();
            return string.Join (string.Empty, rcPattern.Select (Patterns.Complement));
        }
    }
}
