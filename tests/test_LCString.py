import unittest
from goombay import LongestCommonSubstring
from goombay import LongestCommonSubstringMSA


class TestLCSub(unittest.TestCase):
    """Test suite for Longest Common Substring msa_algorithm"""

    def setUp(self):
        """Initialize msa_algorithm for tests"""
        self.algorithm = LongestCommonSubstring()
        self.msa_algorithm = LongestCommonSubstringMSA()

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        # Pairwise
        self.assertEqual(self.algorithm.align("ACTG", "ACTG"), ["ACTG"])
        self.assertEqual(self.algorithm.similarity("ACTG", "ACTG"), 4)
        self.assertEqual(self.algorithm.distance("ACTG", "ACTG"), 0)
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)

        # MSA
        self.assertEqual(self.msa_algorithm.align(["ACTG", "ACTG"]), ["ACTG"])
        self.assertEqual(self.msa_algorithm.similarity(["ACTG", "ACTG"]), 4)
        self.assertEqual(self.msa_algorithm.distance(["ACTG", "ACTG"]), 0)
        self.assertEqual(
            self.msa_algorithm.normalized_similarity(["ACTG", "ACTG"]), 1.0
        )
        self.assertEqual(self.msa_algorithm.normalized_distance(["ACTG", "ACTG"]), 0.0)

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        # Pairwise
        self.assertEqual(self.algorithm.align("AAAA", "TTTT"), [""])
        self.assertEqual(self.algorithm.similarity("AAAA", "TTTT"), 0)
        self.assertEqual(self.algorithm.distance("AAAA", "TTTT"), 4)
        self.assertEqual(self.algorithm.normalized_similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.normalized_distance("AAAA", "TTTT"), 1.0)

        # MSA
        self.assertEqual(self.msa_algorithm.align(["AAAA", "TTTT"]), [""])
        self.assertEqual(self.msa_algorithm.similarity(["AAAA", "TTTT"]), 0)
        self.assertEqual(self.msa_algorithm.distance(["AAAA", "TTTT"]), 4)
        self.assertEqual(
            self.msa_algorithm.normalized_similarity(["AAAA", "TTTT"]), 0.0
        )
        self.assertEqual(self.msa_algorithm.normalized_distance(["AAAA", "TTTT"]), 1.0)

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            (["", "ACTG"], 4),
            (["ACTG", ""], 4),
            (["", ""], 0),
        ]
        for seqs, dist in test_cases:
            with self.subTest(seqs=seqs):
                # Pairwise
                self.assertEqual(self.algorithm.align(seqs[0], seqs[1]), [""])
                self.assertEqual(self.algorithm.similarity(seqs[0], seqs[1]), 0)
                self.assertEqual(self.algorithm.distance(seqs[0], seqs[1]), dist)

                # MSA
                self.assertEqual(self.msa_algorithm.align(seqs), [""])
                self.assertEqual(self.msa_algorithm.similarity(seqs), 0)
                self.assertEqual(self.msa_algorithm.distance(seqs), dist)

    def test_single_character_sequences(self):
        """Test behavior with single character sequences"""
        # Pairwise
        self.assertEqual(self.algorithm.align("A", "A"), [""])
        self.assertEqual(self.algorithm.similarity("A", "A"), 0)
        self.assertEqual(self.algorithm.normalized_similarity("A", "A"), 1)
        self.assertEqual(self.algorithm.distance("A", "A"), 0)
        self.assertEqual(self.algorithm.normalized_distance("A", "A"), 0)

        self.assertEqual(self.algorithm.align("A", "T"), [""])
        self.assertEqual(self.algorithm.similarity("A", "T"), 0)
        self.assertEqual(self.algorithm.normalized_similarity("A", "T"), 0)
        self.assertEqual(self.algorithm.distance("A", "T"), 1)
        self.assertEqual(self.algorithm.normalized_distance("A", "T"), 1)

        # MSA
        self.assertEqual(self.msa_algorithm.align(["A", "A"]), [""])
        self.assertEqual(self.msa_algorithm.similarity(["A", "A"]), 0)
        self.assertEqual(self.msa_algorithm.normalized_similarity(["A", "A"]), 1)
        self.assertEqual(self.msa_algorithm.distance(["A", "A"]), 0)
        self.assertEqual(self.msa_algorithm.normalized_distance(["A", "A"]), 0)

        self.assertEqual(self.msa_algorithm.align(["A", "T"]), [""])
        self.assertEqual(self.msa_algorithm.similarity(["A", "T"]), 0)
        self.assertEqual(self.msa_algorithm.normalized_similarity(["A", "T"]), 0)
        self.assertEqual(self.msa_algorithm.distance(["A", "T"]), 1)
        self.assertEqual(self.msa_algorithm.normalized_distance(["A", "T"]), 1)

    def test_case_sensitivity(self):
        """Test that matching is case-sensitive (unlike LCS)"""
        # Pairwise
        self.assertEqual(self.algorithm.align("ACTG", "actg"), ["ACTG"])
        self.assertEqual(self.algorithm.similarity("ACTG", "actg"), 4)
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "actg"), 1)
        self.assertEqual(self.algorithm.distance("ACTG", "actg"), 0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "actg"), 0)

        # MSA
        self.assertEqual(self.msa_algorithm.align(["ACTG", "actg"]), ["ACTG"])
        self.assertEqual(self.msa_algorithm.similarity(["ACTG", "actg"]), 4)
        self.assertEqual(self.msa_algorithm.normalized_similarity(["ACTG", "actg"]), 1)
        self.assertEqual(self.msa_algorithm.distance(["ACTG", "actg"]), 0)
        self.assertEqual(self.msa_algorithm.normalized_distance(["ACTG", "actg"]), 0)

    def test_longest_shared_substring(self):
        """Test known shared substrings"""
        # Pairwise
        self.assertEqual(self.algorithm.align("GATTACA", "TTAC"), ["TTAC"])
        self.assertEqual(self.algorithm.similarity("GATTACA", "TTAC"), 4)
        self.assertEqual(self.algorithm.normalized_similarity("GATTACA", "TTAC"), 1)
        self.assertEqual(self.algorithm.normalized_distance("GATTACA", "TTAC"), 0)

        # MSA
        self.assertEqual(self.msa_algorithm.align(["GATTACA", "TTAC"]), ["TTAC"])
        self.assertEqual(self.msa_algorithm.similarity(["GATTACA", "TTAC"]), 4)
        self.assertEqual(
            self.msa_algorithm.normalized_similarity(["GATTACA", "TTAC"]), 1
        )
        self.assertEqual(self.msa_algorithm.normalized_distance(["GATTACA", "TTAC"]), 0)

    def test_different_lengths(self):
        """Test sequences of different lengths"""
        test_cases = [
            (["ABCDE", "XABCY"], ["ABC"]),
            (["HELLO", "YELLOWS"], ["ELLO"]),
            (["STRAWBERRY", "RAWR"], ["RAW"]),
            (["HUMAN", "CHIMPANZEE"], ["AN"]),
        ]
        for seqs, expected in test_cases:
            with self.subTest(seqs=seqs):
                # Pairwise
                self.assertEqual(expected, self.algorithm.align(seqs[0], seqs[1]))
                self.assertEqual(
                    self.algorithm.similarity(seqs[0], seqs[1]), len(expected[0])
                )

                # MSA
                self.assertEqual(expected, self.msa_algorithm.align(seqs))
                self.assertEqual(self.msa_algorithm.similarity(seqs), len(expected[0]))

    def test_multiple_sequences(self):
        """Test alignment with more than 2 sequences"""
        sequences = ["CCTAGGAC", "GGACCTAG", "TAGGAC"]
        result = self.msa_algorithm.align(sequences)
        self.assertTrue(all(sub in seq for seq in sequences for sub in result))
        self.assertTrue(all(len(s) == len(result[0]) for s in result))

    def test_substring_vs_subsequence(self):
        """Distinguish from LCS: only substrings, no gaps allowed"""
        # Pairwise
        self.assertEqual(self.algorithm.align("ABCDEF", "ACE"), [""])
        self.assertEqual(self.algorithm.similarity("ABCDEF", "ACE"), 0)

        # MSA
        self.assertEqual(self.msa_algorithm.align(["ABCDEF", "ACE"]), [""])
        self.assertEqual(self.msa_algorithm.similarity(["ABCDEF", "ACE"]), 0)

    def test_similarity_distance_consistency(self):
        """Check that similarity and distance always add up correctly"""
        test_cases = [
            ["ACTG", "CTGA"],
            ["GATTACA", "TTACAGG"],
            ["AAAA", "TTTT"],
        ]
        for seqs in test_cases:
            with self.subTest(seqs=seqs):
                max_len = len(max(seqs, key=len))
                sim = self.msa_algorithm.similarity(seqs)
                dist = self.msa_algorithm.distance(seqs)
                self.assertEqual(dist, max_len - sim)

    def test_multiple_sequences_shared_substring(self):
        """Test longest common substring across multiple sequences"""
        sequences = ["GATTACA", "TTACAGG", "CTTACAA", "ACTTACATTGCA"]
        result = self.msa_algorithm.align(sequences)
        # The substring "TTACA" is shared across all
        self.assertEqual(["TTACA"], result)
        self.assertTrue(all("TTACA" in seq for seq in sequences))
        self.assertEqual(self.msa_algorithm.similarity(sequences), 5)
        self.assertEqual(
            self.msa_algorithm.distance(sequences), len(max(sequences, key=len)) - 5
        )


if __name__ == "__main__":
    unittest.main()
