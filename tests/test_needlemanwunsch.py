import unittest
import numpy
from goombay import NeedlemanWunsch

class TestNeedlemanWunsch(unittest.TestCase):
    """Test suite for Needleman-Wunsch global alignment algorithm"""

    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = NeedlemanWunsch()

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        # Test alignment
        self.assertEqual(self.algorithm.align("ACTG", "ACTG"), "ACTG\nACTG")

        # Test scoring
        self.assertEqual(self.algorithm.similarity("ACTG", "ACTG"), 4 * self.algorithm.match_score)
        self.assertEqual(self.algorithm.distance("ACTG", "ACTG"), 0.0)

        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        # Test alignment
        self.assertEqual(self.algorithm.align("AAAA", "TTTT"), "AAAA\nTTTT")

        # Test scoring
        self.assertEqual(self.algorithm.similarity("AAAA", "TTTT"), -4 * self.algorithm.mismatch_penalty)
        self.assertEqual(self.algorithm.distance("AAAA", "TTTT"), 4 * self.algorithm.mismatch_penalty)

        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.normalized_distance("AAAA", "TTTT"), 1.0)

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "ACTG"),    # Empty query
            ("ACTG", ""),    # Empty subject
            ("", "")         # Both empty
        ]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                if query == subject == "":
                    self.assertEqual(self.algorithm.align(query, subject), "\n")
                    self.assertEqual(self.algorithm.similarity(query, subject), 1.0)
                    self.assertEqual(self.algorithm.distance(query, subject), 0.0)
                else:
                    non_empty = query or subject
                    expected_gaps = "-" * len(non_empty)
                    self.assertEqual(
                        self.algorithm.align(query, subject),
                        f"{expected_gaps if not query else query}\n{expected_gaps if not subject else subject}"
                    )
                    self.assertEqual(self.algorithm.similarity(query, subject), -2.0 * len(non_empty))
                    self.assertEqual(self.algorithm.distance(query, subject), 2.0 * len(non_empty))

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.align("A", "A"), "A\nA")
        self.assertEqual(self.algorithm.similarity("A", "A"), self.algorithm.match_score)
        self.assertEqual(self.algorithm.distance("A", "A"), 0.0)

        # Test mismatch
        self.assertEqual(self.algorithm.align("A", "T"), "A\nT")
        self.assertEqual(self.algorithm.similarity("A", "T"), -self.algorithm.mismatch_penalty)
        self.assertEqual(self.algorithm.distance("A", "T"), self.algorithm.mismatch_penalty)

    def test_case_sensitivity(self):
        """Test that matching is case-insensitive"""
        test_cases = [
            ("ACTG", "actg"),
            ("AcTg", "aCtG"),
            ("actg", "ACTG")
        ]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.align(query, subject),
                    self.algorithm.align(query.upper(), subject.upper())
                )
                self.assertEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(query.upper(), subject.upper())
                )

    def test_matrix_shape(self):
        """Test matrix dimensions"""
        query = "ACTG"
        subject = "ACT"
        score, pointer = self.algorithm(query, subject)
        expected_shape = (len(query) + 1, len(subject) + 1)
        self.assertEqual(score.shape, expected_shape)
        self.assertEqual(pointer.shape, expected_shape)

    def test_matrix_values(self):
        """Test specific matrix values"""
        score, pointer = self.algorithm("AC", "AT")

        expected_score = numpy.array([
            [ 0, -2, -4],
            [-2,  2,  0],
            [-4,  0,  1]
        ])
        numpy.testing.assert_array_equal(score, expected_score)

        # Test pointer values are valid
        self.assertTrue(numpy.all(pointer >= 0))
        self.assertTrue(numpy.all(pointer <= 7))  # Max valid pointer value

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTG", "ACT", "ACTG\nACT-"),    # Longer query
            ("ACT", "ACTG", "ACT-\nACTG"),    # Longer subject
            ("ACGT", "AGT", "ACGT\nA-GT")     # Internal gap
        ]

        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

    def test_scoring_parameters(self):
        """Test behavior with different scoring parameters"""
        custom_algorithm = NeedlemanWunsch(
            match_score=1,
            mismatch_penalty=2,
            gap_penalty=3
        )

        # Test alignment
        self.assertEqual(custom_algorithm.align("ACGT", "AGT"), "ACGT\nA-GT")

        # Test matrix values
        score, _ = custom_algorithm("AC", "AT")
        expected_score = numpy.array([
            [ 0, -3, -6],
            [-3,  1, -2],
            [-6, -2, -1]
        ])
        numpy.testing.assert_array_equal(score, expected_score)

    def test_normalization(self):
        """Test normalization behavior"""
        test_cases = [
            ("ACTG", "ACTG", 1.0, 0.0),    # identical
            ("ACTG", "AATG", 0.75, 0.25),  # one mismatch
            ("ACTG", "AAAG", 0.5, 0.5),    # two mismatches
            ("ACTG", "AAAA", 0.25, 0.75),  # three mismatches
            ("ACTG", "VVVV", 0.0, 1.0)     # all mismatches
        ]

        for query, subject, exp_sim, exp_dist in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(self.algorithm.normalized_similarity(query, subject), exp_sim)
                self.assertAlmostEqual(self.algorithm.normalized_distance(query, subject), exp_dist)

if __name__ == '__main__':
    unittest.main()
