import unittest
import numpy
from goombay import LowranceWagner


class TestLowranceWagner(unittest.TestCase):
    """Test suite for Lowrance-Wagner edit distance algorithm"""

    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = LowranceWagner()

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "", "\n"),  # Both empty
            ("", "ACTG", "----\nACTG"),  # Empty query
            ("ACTG", "", "ACTG\n----"),  # Empty subject
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

    def test_single_character(self):
        """Test behavior with single character sequences"""
        test_cases = [
            ("A", "A", "A\nA"),  # Same character
            ("A", "T", "A\nT"),  # Different characters
            ("T", "A", "T\nA"),  # Different characters, reversed
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

    def test_dna_sequences(self):
        """Test behavior with DNA sequences"""
        test_cases = [
            ("ACTG", "ACTG", "ACTG\nACTG"),  # Identical
            ("AAAA", "TTTT", "AAAA\nTTTT"),  # All different
            ("AACTG", "ACTGG", "AACT-G\n-ACTGG"),  # Nested overlap
            ("AGCT", "TAGC", "-AGCT\nTAGC-"),  # Partial overlap
            ("CCTT", "GGGA", "CCTT\nGGGA"),  # No overlap
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

    def test_transposition_alignments(self):
        """Test alignments involving transpositions"""
        test_cases = [
            ("ABCD", "BADC", "ABCD\nBADC"),  # Two clear transpositions
            ("ACGT", "CATG", "ACGT\nCATG"),  # Two clear transpositions
            ("ATCG", "TACG", "ATCG\nTACG"),  # One transposition
            ("AABB", "BBAA", "AABB\nBBAA"),  # One transposition
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

    def test_transposition_distances(self):
        """Test distances involving transpositions"""
        test_cases = [
            ("ABCD", "BADC", 2),  # Two transpositions
            ("ACGT", "CATG", 2),  # Two transpositions
            ("ATCG", "TACG", 1),  # One transposition
            ("AABB", "BBAA", 3),  # One transposition + two substitutions
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.distance(query, subject), expected)

    def test_similarity(self):
        """Test similarity calculation"""
        test_cases = [
            ("ACTG", "ACTG", 4.0),  # Identical - 4 matches
            ("AAAA", "TTTT", 0.0),  # All mismatches
            ("AACTG", "ACTGG", 3.0),  # 3 matches, 2 gaps
            ("", "", 1.0),  # Empty sequences
            ("ACTG", "", 0.0),  # One empty sequence
            ("A", "T", 0.0),  # One mismatch
            ("GGCC", "CCGG", 1.0),  # One transposition
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.similarity(query, subject), expected)

    def test_distance(self):
        """Test distance calculation"""
        test_cases = [
            ("ACTG", "ACTG", 0),  # Identical
            ("AAAA", "TTTT", 4),  # All different
            ("AACTG", "ACTGG", 2),  # One insertion, one deletion
            ("", "", 0),  # Empty sequences
            ("ACTG", "", 4),  # One empty sequence
            ("A", "T", 1),  # Single substitution
            ("AATCG", "ATCGG", 2),  # Two operations needed
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.distance(query, subject), expected)

    def test_matrix(self):
        """Test matrix calculation"""
        test_cases = [
            # Empty sequences
            ("", "", numpy.zeros((1, 1))),
            # Single character match
            ("A", "A", numpy.array([[0, 1], [1, 0]])),
            # Single character mismatch
            ("A", "T", numpy.array([[0, 1], [1, 1]])),
            # Simple sequence
            ("AT", "TA", numpy.array([[0, 1, 2], [1, 1, 1], [2, 1, 1]])),
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                result = self.algorithm.matrix(query, subject)
                numpy.testing.assert_array_equal(result, expected)

    def test_normalized_similarity(self):
        """Test normalized similarity calculation"""
        test_cases = [
            ("ACTG", "ACTG", 1.0),  # Identical
            ("AAAA", "TTTT", 0.0),  # All different
            ("AACTG", "ACTGG", 0.6),  # 3 matches out of 5
            ("", "", 1.0),  # Empty sequences
            ("ACTG", "", 0.0),  # One empty sequence
            ("A", "T", 0.0),  # Single mismatch
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.normalized_similarity(query, subject),
                    expected,
                    places=3,
                )

    def test_normalized_distance(self):
        """Test normalized distance calculation"""
        test_cases = [
            ("ACTG", "ACTG", 0.0),  # Identical
            ("AAAA", "TTTT", 1.0),  # All different
            ("AACTG", "ACTGG", 0.4),  # 2 operations out of 5
            ("", "", 0.0),  # Empty sequences
            ("ACTG", "", 1.0),  # One empty sequence
            ("A", "T", 1.0),  # Single mismatch
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.normalized_distance(query, subject),
                    expected,
                    places=3,
                )

    def test_all_alignments(self):
        """Test returning multiple optimal alignments"""
        test_cases = [
            ("ACCG", "ACG", ["ACCG\nAC-G", "ACCG\nA-CG"], 2),
            (
                "ATGTGTA",
                "ATA",
                ["ATGTGTA\nAT----A", "ATGTGTA\nA----TA", "ATGTGTA\nA--T--A"],
                3,
            ),
            ("ACGGCT", "ACT", ["ACGGCT\nAC---T", "ACGGCT\nA---CT"], 2),
            ("CCGA", "CG", ["CCGA\nC-G-", "CCGA\n-CG-"], 2),
        ]
        for query, subject, alignments, length in test_cases:
            with self.subTest(query=query, subject=subject):
                res = self.algorithm.align(query, subject, all_alignments=True)
                self.assertEqual(length, len(res))
                for alignment in alignments:
                    self.assertIn(alignment, res)


if __name__ == "__main__":
    unittest.main()
