import unittest
import numpy
from goombay import SmithWaterman


class TestSmithWaterman(unittest.TestCase):
    """Test suite for Smith-Waterman local alignment algorithm"""

    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = SmithWaterman()

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        # Test alignment
        self.assertEqual(self.algorithm.align("ACTG", "ACTG"), "ACTG\nACTG")

        # Test scoring
        self.assertEqual(
            self.algorithm.similarity("ACTG", "ACTG"), 4 * self.algorithm.match_score
        )
        self.assertEqual(self.algorithm.distance("ACTG", "ACTG"), 0)

        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        # Test alignment
        self.assertEqual(
            self.algorithm.align("AAAA", "TTTT"), "There is no local alignment!"
        )

        # Test scoring
        self.assertEqual(self.algorithm.similarity("AAAA", "TTTT"), 0.0)

        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.normalized_distance("AAAA", "TTTT"), 1.0)

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "ACTG", 4, 1, 0, 0),  # Empty query
            ("ACTG", "", 4, 1, 0, 0),  # Empty subject
            ("", "", 0, 0, 1, 1),  # Empty subject
        ]

        for query, subject, dist, norm_dist, sim, norm_sim in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.align(query, subject), "There is no local alignment!"
                )
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(
                    self.algorithm.normalized_similarity(query, subject), norm_sim
                )
                self.assertEqual(self.algorithm.distance(query, subject), dist)
                self.assertEqual(
                    self.algorithm.normalized_distance(query, subject), norm_dist
                )

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.align("A", "A"), "A\nA")
        self.assertEqual(self.algorithm.similarity("A", "A"), 1.0)
        self.assertEqual(self.algorithm.distance("A", "A"), 0.0)

        # Test mismatch
        self.assertEqual(self.algorithm.align("A", "T"), "There is no local alignment!")
        self.assertEqual(self.algorithm.similarity("A", "T"), 0.0)
        self.assertEqual(self.algorithm.distance("A", "T"), 1.0)

    def test_case_sensitivity(self):
        """Test that matching is case-insensitive"""
        test_cases = [("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.align(query, subject),
                    self.algorithm.align(query.upper(), subject.upper()),
                )
                self.assertEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(query.upper(), subject.upper()),
                )

    def test_matrix_shape(self):
        """Test matrix dimensions"""
        query = "ACTG"
        subject = "ACT"
        score = self.algorithm(query, subject)
        expected_shape = (len(query) + 1, len(subject) + 1)
        self.assertEqual(score.shape, expected_shape)

    def test_matrix_values(self):
        """Test specific matrix values"""
        score = self.algorithm("AC", "AT")

        expected_score = numpy.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
        numpy.testing.assert_array_equal(score, expected_score)

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTGACTG", "ACT", "ACT\nACT", 3, 5),  # Longer query
            ("ACT", "ACTGACTG", "ACT\nACT", 3, 5),  # Longer subject
            ("ACGTACGT", "AGT", "GT\nGT", 2, 6),  # Internal match
        ]

        for query, subject, expected, sim, dist in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(self.algorithm.distance(query, subject), dist)

    def test_scoring_parameters(self):
        """Test behavior with different scoring parameters"""
        custom_algorithm = SmithWaterman(
            match_score=1, mismatch_penalty=2, gap_penalty=3
        )

        # Test alignment
        self.assertEqual(custom_algorithm.align("ACGT", "AGT"), "GT\nGT")

        # Test matrix values
        score = custom_algorithm("AC", "AT")
        expected_score = numpy.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
        numpy.testing.assert_array_equal(score, expected_score)

    def test_local_alignment_behavior(self):
        """Test specific local alignment behaviors"""
        test_cases = [
            ("CGATC", "GTATG", "AT\nAT", 2, 3),  # Find internal match
            ("AAAGGGGTTT", "AAATTT", "TTT\nTTT", 3, 7),  # Multiple possible alignments
            ("ACGTACGT", "TACGTAC", "ACGTAC\nACGTAC", 6, 2),  # Longer local alignment
        ]

        for query, subject, expected, sim, dist in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(self.algorithm.distance(query, subject), dist)


if __name__ == "__main__":
    unittest.main()
