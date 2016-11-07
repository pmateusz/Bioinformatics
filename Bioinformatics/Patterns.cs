using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Bioinformatics
{
    class Patterns
    {
        public static readonly string[] Nucleotides = new[] { "A", "C", "T", "G" };

        public static readonly char[] RawNucleotides = new[] { 'A', 'C', 'T', 'G' };

        public static readonly Dictionary<char, int> Peptides = new Dictionary<char, int> () {
            {'G',57},
            {'A',71},
            {'S',87},
            {'P',97},
            {'V',99},
            {'T',101},
            {'C',103},
            {'I',113},
            {'L',113},
            {'N',114},
            {'D',115},
            {'K',128},
            {'Q',128},
            {'E',129},
            {'M',131},
            {'H',137},
            {'F',147},
            {'R',156},
            {'Y',163},
            {'W',186}
        };

        public static readonly Dictionary<string, char> Codons = new Dictionary<string, char> () {
            {"AAA",'K'},
            {"AAC",'N'},
            {"AAG",'K'},
            {"AAU",'N'},
            {"ACA",'T'},
            {"ACC",'T'},
            {"ACG",'T'},
            {"ACU",'T'},
            {"AGA",'R'},
            {"AGC",'S'},
            {"AGG",'R'},
            {"AGU",'S'},
            {"AUA",'I'},
            {"AUC",'I'},
            {"AUG",'M'},
            {"AUU",'I'},
            {"CAA",'Q'},
            {"CAC",'H'},
            {"CAG",'Q'},
            {"CAU",'H'},
            {"CCA",'P'},
            {"CCC",'P'},
            {"CCG",'P'},
            {"CCU",'P'},
            {"CGA",'R'},
            {"CGC",'R'},
            {"CGG",'R'},
            {"CGU",'R'},
            {"CUA",'L'},
            {"CUC",'L'},
            {"CUG",'L'},
            {"CUU",'L'},
            {"GAA",'E'},
            {"GAC",'D'},
            {"GAG",'E'},
            {"GAU",'D'},
            {"GCA",'A'},
            {"GCC",'A'},
            {"GCG",'A'},
            {"GCU",'A'},
            {"GGA",'G'},
            {"GGC",'G'},
            {"GGG",'G'},
            {"GGU",'G'},
            {"GUA",'V'},
            {"GUC",'V'},
            {"GUG",'V'},
            {"GUU",'V'},
            {"UAA",'?'},
            {"UAC",'Y'},
            {"UAG",'?'},
            {"UAU",'Y'},
            {"UCA",'S'},
            {"UCC",'S'},
            {"UCG",'S'},
            {"UCU",'S'},
            {"UGA",'?'},
            {"UGC",'C'},
            {"UGG",'W'},
            {"UGU",'C'},
            {"UUA",'L'},
            {"UUC",'F'},
            {"UUG",'L'},
            {"UUU",'F'}};

        public static readonly Dictionary<char, string[]> RawCodons = new Dictionary<char, string[]> () {
            {'A', new []{"GCA","GCC","GCG","GCU"}},
            {'C', new []{"UGC","UGU"}},
            {'D', new []{"GAC","GAU"}},
            {'E', new []{"GAA","GAG"}},
            {'F', new []{"UUC","UUU"}},
            {'G', new []{"GGA","GGC","GGG","GGU"}},
            {'H', new []{"CAC","CAU"}},
            {'I', new []{"AUA","AUC","AUU"}},
            {'K', new []{"AAA","AAG"}},
            {'L', new []{"CUA","CUU","CUC","CUG","UUA","UUG"}},
            {'M', new []{"AUG"}},
            {'N', new []{"AAC","AAU"}},
            {'P', new []{"CCA","CCC","CCG","CCU"}},
            {'Q', new []{"CAA","CAG"}},
            {'R', new []{"AGA","AGG","CGA","CGC","CGG","CGU"}},
            {'S', new []{"AGC","AGU","UCA","UCC","UCG","UCU"}},
            {'T', new []{"ACA","ACC","ACG","ACU"}},
            {'V', new []{"GUA", "GUC", "GUG", "GUU"}},
            {'Y', new []{"UAC","UAU"}},
            {'?', new []{"UAG","UAA","UGA"}},
            {'W', new []{"UGG"}}
        };

        public static string TranslateNucleotides (string code)
        {
            const int Length = 3;

            var translation = new StringBuilder ();
            for (var position = 0; position < code.Length; position += Length)
            {
                var codon = code.Substring (position, Length);
                translation.Append (Codons[codon]);
            }
            return translation.ToString ();
        }

        public static char Nucleotide (string codon)
        {
            return Codons[codon];
        }

        public static char Complement (char nucleotide)
        {
            switch (nucleotide)
            {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
            }

            throw new ArgumentException ();
        }

        public static string NumberToPattern (int number, int k)
        {
            var currentNumber = number;
            var pattern = new List<char> ();
            while (k > 0)
            {
                pattern.Add (ToSymbol (currentNumber % 4));
                currentNumber = currentNumber / 4;
                --k;
            }
            pattern.Reverse ();
            return string.Join (string.Empty, pattern);
        }

        public static IEnumerable<char> FilterNucleotides (char except)
        {
            foreach (var nucleotide in Patterns.RawNucleotides)
            {
                if (nucleotide == except)
                {
                    continue;
                }

                yield return nucleotide;
            }
        }

        public static int ToNumber (string pattern)
        {
            if (pattern.Length == 0)
            {
                return 0;
            }

            var result = 0;
            foreach (var symbol in pattern)
            {
                result = 4 * result + ToNumber (symbol);
            }
            return result;
        }

        public static int ToNumber (char symbol)
        {
            switch (symbol)
            {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
            }

            throw new ArgumentException ();
        }

        public static char ToSymbol (int number)
        {
            switch (number)
            {
                case 0: return 'A';
                case 1: return 'C';
                case 2: return 'G';
                case 3: return 'T';
            }

            throw new ArgumentException ();
        }
    }
}
