import unittest
import numpy
from goombay import WatermanSmithBeyer


class TestWatermanSmithBeyer(unittest.TestCase):
    """Test suite for Waterman-Smith-Beyer alignment algorithm"""

    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = WatermanSmithBeyer()

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        score, pointer = self.algorithm("ACTG", "ACTG")
        # Test scoring
        self.assertEqual(score[-1][-1], 4 * self.algorithm.match)
        self.assertEqual(pointer[-1][-1], (2, 0, 0))  # All diagonal moves

        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        score, pointer = self.algorithm("AAAA", "TTTT")
        # Test scoring
        self.assertEqual(score[-1][-1], -4 * self.algorithm.mismatch)  # All mismatches

        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.normalized_distance("AAAA", "TTTT"), 1.0)

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "ACTG"),  # Empty query
            ("ACTG", ""),  # Empty subject
            ("", ""),  # Both empty
        ]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                score, _ = self.algorithm(query, subject)
                expected_score = -self.algorithm.gap - self.algorithm.continued_gap * (
                    len(query or subject)
                )  # gap + continue_gaps
                if query == subject == "":
                    self.assertEqual(self.algorithm.similarity(query, subject), 1)
                elif not subject:
                    self.assertEqual(score[len(query)][0], expected_score)
                elif not query:
                    self.assertEqual(score[0][len(subject)], expected_score)

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        score, pointer = self.algorithm("A", "A")
        self.assertEqual(score[-1][-1], self.algorithm.match)
        self.assertEqual(pointer[-1][-1], (2, 0, 0))  # Diagonal move

        # Test mismatch
        score, pointer = self.algorithm("A", "T")
        self.assertEqual(score[-1][-1], -self.algorithm.mismatch)
        self.assertEqual(pointer[-1][-1], (2, 0, 0))  # Diagonal move

    def test_case_sensitivity(self):
        """Test that matching is case-insensitive"""
        test_cases = [("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                score1, _ = self.algorithm(query, subject)
                score2, _ = self.algorithm(query.upper(), subject.upper())
                numpy.testing.assert_array_equal(score1, score2)

    def test_matrix_shape(self):
        """Test matrix dimensions"""
        query, subject = "ACTG", "AGT"
        score, pointer = self.algorithm(query, subject)
        expected_shape = (len(query) + 1, len(subject) + 1)
        self.assertEqual(score.shape, expected_shape)
        self.assertEqual(pointer.shape, expected_shape)

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTG", "ACT"),  # Longer query
            ("ACT", "ACTG"),  # Longer subject
            ("ACGT", "AGT"),  # Internal gap
            ("HOLYROMANEMPIRE", "HOLYPIRE"),  # Longer subject
        ]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                score, pointer = self.algorithm(query, subject)
                self.assertEqual(score.shape, (len(query) + 1, len(subject) + 1))

    def test_gap_penalties(self):
        """Test behavior with gaps"""
        test_cases = [
            ("AT", "A-T"),  # Single gap
            ("ACTG", "AG"),  # Multiple gaps
            ("A", "A--A"),  # Consecutive gaps
        ]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                _, pointer_matrix = self.algorithm(query, subject)
                has_gap = any(
                    pointer[0] in (3, 4) for row in pointer_matrix for pointer in row
                )
                self.assertTrue(has_gap)

    def test_scoring_parameters(self):
        """Test behavior with different scoring parameters"""
        custom_algorithm = WatermanSmithBeyer(
            match=1, mismatch=2, new_gap=3, continued_gap=1
        )

        score, _ = custom_algorithm("AC", "AT")
        expected_score = numpy.array(
            [[0.0, -4.0, -5.0], [-4.0, 1.0, -3.0], [-5.0, -3.0, -1.0]]
        )
        numpy.testing.assert_array_almost_equal(score, expected_score)

    def test_all_alignments(self):
        """Test returning multiple optimal alignments"""
        test_cases = [
            ("ACCG", "ACG", ["ACCG\nAC-G", "ACCG\nA-CG"], 2),
            (
                "ATGTGTA",
                "ATA",
                ["ATGTGTA\nAT----A", "ATGTGTA\nA----TA"],
                2,
            ),
            ("ACGGCT", "ACT", ["ACGGCT\nAC---T", "ACGGCT\nA---CT"], 2),
            ("CCGA", "CG", ["CCGA\nCG--", "CCGA\nC--G"], 2),
        ]
        for query, subject, alignments, length in test_cases:
            with self.subTest(query=query, subject=subject):
                res = self.algorithm.align(query, subject, all_alignments=True)
                self.assertEqual(length, len(res))
                for alignment in alignments:
                    self.assertIn(alignment, res)

    def test_gap_functions_affine(self):
        """Test affine gap function produces expected alignment"""
        affine = WatermanSmithBeyer(
            match=2, mismatch=1, new_gap=3, continued_gap=1, gap_function="affine"
        )
        query, subject = "ACCGT", "CT"
        alignment = affine.align(query, subject)
        # Affine gap: expect one gap, not two single gaps
        self.assertIn(alignment, ["ACCGT\n---CT", "ACCGT\nC---T"])

    def test_gap_functions_quadratic(self):
        """Test quadratic gap function penalizes longer gaps more heavily"""
        quad = WatermanSmithBeyer(
            match=2, mismatch=1, new_gap=3, continued_gap=1, gap_function="quadratic"
        )
        query, subject = "ACCGT", "GT"
        alignment = quad.align(query, subject)
        self.assertIn(alignment, ["ACCGT\n---GT"])

    def test_gap_functions_logarithmic(self):
        """Test logarithmic gap function penalizes longer gaps less than quadratic"""
        log = WatermanSmithBeyer(
            match=2, mismatch=1, new_gap=3, continued_gap=1, gap_function="logarithmic"
        )
        query, subject = "ACCGT", "CT"
        alignment = log.align(query, subject)
        self.assertIn(alignment, ["ACCGT\n---CT", "ACCGT\nC---T"])

    def test_gap_function_invalid(self):
        """Test that invalid gap function raises ValueError"""
        with self.assertRaises(ValueError):
            self.algorithm.gap_function = "invalid"
            self.algorithm("ACGT", "AGT")


if __name__ == "__main__":
    unittest.main()
