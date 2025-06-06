import unittest
import numpy
from goombay import Jaro


class TestJaro(unittest.TestCase):
    """Test suite for Jaro distance algorithm"""

    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = Jaro()

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        test_cases = [
            "ACTG",  # DNA sequence
            "1010",  # Binary string
            "Hello",  # Regular string
            "A",  # Single character
            "",  # Empty string
        ]

        for sequence in test_cases:
            with self.subTest(sequence=sequence):
                # Test similarity
                self.assertEqual(self.algorithm.similarity(sequence, sequence), 1.0)

                # Test distance
                self.assertEqual(self.algorithm.distance(sequence, sequence), 0.0)

                # Test normalization
                self.assertEqual(
                    self.algorithm.normalized_similarity(sequence, sequence), 1.0
                )
                self.assertEqual(
                    self.algorithm.normalized_distance(sequence, sequence), 0.0
                )

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        test_cases = [
            ("AAAA", "TTTT"),  # DNA sequences
            ("0000", "1111"),  # Binary strings
            ("Hello", "Warvd"),  # Regular strings
            ("ABC", "XYZ"),  # No matching characters
        ]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                # Test similarity
                self.assertEqual(self.algorithm.similarity(query, subject), 0.0)

                # Test distance
                self.assertEqual(self.algorithm.distance(query, subject), 1.0)

                # Test normalization
                self.assertEqual(
                    self.algorithm.normalized_similarity(query, subject), 0.0
                )
                self.assertEqual(
                    self.algorithm.normalized_distance(query, subject), 1.0
                )

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "", 0, 1, "\n"),  # Both empty sequences
            ("A", "", 1, 0, "A\n-"),  # One empty subject sequence
            ("", "A", 1, 0, "-\nA"),  # One empty query sequence
            ("", "ACTG", 1, 0, "----\nACTG"),  # Longer empty query sequence
            ("ACTG", "", 1, 0, "ACTG\n----"),  # Longer empty subject sequence
        ]

        for query, subject, dist, sim, aligned in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(self.algorithm.distance(query, subject), dist)
                self.assertEqual(self.algorithm.align(query, subject), aligned)

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.similarity("A", "A"), 1.0)
        self.assertEqual(self.algorithm.distance("A", "A"), 0.0)
        self.assertEqual(self.algorithm.align("A", "A"), "A\nA")

        # Test mismatch
        self.assertEqual(self.algorithm.similarity("A", "T"), 0.0)
        self.assertEqual(self.algorithm.distance("A", "T"), 1.0)
        self.assertEqual(self.algorithm.align("A", "T"), "A-\n-T")

    def test_case_sensitivity(self):
        """Test that matching is case-insensitive"""
        test_cases = [("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.similarity(query, subject), 1.0)
                self.assertEqual(self.algorithm.distance(query, subject), 0.0)

    def test_transpositions(self):
        """Test handling of character transpositions"""
        test_cases = [
            ("MARTHA", "MARHTA", 0.944),  # One transposition
            ("DIXON", "DICKSONX", 0.767),  # Multiple differences
            ("JELLYFISH", "SMELLYFISH", 0.896),  # Prefix difference
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.similarity(query, subject), expected, places=3
                )

    def test_matrix(self):
        """Test matrix output for matches"""
        test_cases = [
            (
                "CRATE",
                "TRACE",
                numpy.array(
                    [
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1, 0],
                        [0, 1, 1, 1, 1, 2],
                        [0, 1, 2, 2, 2, 2],
                        [0, 1, 1, 2, 2, 2],
                        [0, 1, 2, 2, 3, 3],
                    ]
                ),
            ),
            (
                "MARTHA",
                "MARHTA",
                numpy.array(
                    [
                        [0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 1, 1, 1, 1],
                        [0, 1, 2, 2, 2, 2, 2],
                        [0, 1, 2, 3, 3, 3, 3],
                        [0, 1, 2, 3, 4, 4, 4],
                        [0, 1, 2, 3, 4, 5, 5],
                        [0, 1, 2, 3, 4, 5, 6],
                    ]
                ),
            ),
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                result = self.algorithm.matrix(query, subject)
                numpy.testing.assert_array_equal(result, expected)

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("MARTHA", "MART", 0.889),  # Shorter second string
            ("MART", "MARTHA", 0.889),  # Shorter first string
            ("MARTHA", "MARHTA", 0.944),
            ("DIC", "DICKSON", 0.810),  # Significant length difference
            ("DIXON", "DICKSONX", 0.767),
            ("FAREMVIEL", "FARMVILLE", 0.884),
            ("CRATE", "TRACE", 0.733),
            ("DWAYNE", "DUANE", 0.822),
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.similarity(query, subject), expected, places=3
                )
                self.assertAlmostEqual(
                    self.algorithm.normalized_similarity(query, subject),
                    expected,
                    places=3,
                )


if __name__ == "__main__":
    unittest.main()
