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
        self.assertEqual(score[-1][-1], 4 * self.algorithm.match_score)
        self.assertEqual(pointer[-1][-1], 2)  # All diagonal moves
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        score, pointer = self.algorithm("AAAA", "TTTT")
        # Test scoring
        self.assertEqual(score[-1][-1], -4*self.algorithm.mismatch_penalty)  # All mismatches
        self.assertTrue(numpy.all(pointer[1:,1:] == 2))  # All diagonal moves
        
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
                score, _ = self.algorithm(query, subject)
                expected_score = -4 - (len(query or subject))  # new_gap + continue_gaps
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
        self.assertEqual(score[-1][-1], self.algorithm.match_score)
        self.assertEqual(pointer[-1][-1], 2)  # Diagonal move
        
        # Test mismatch
        score, pointer = self.algorithm("A", "T")
        self.assertEqual(score[-1][-1], -self.algorithm.mismatch_penalty)
        self.assertEqual(pointer[-1][-1], 2)  # Diagonal move

    def test_case_sensitivity(self):
        """Test that matching is case-insensitive"""
        test_cases = [
            ("ACTG", "actg"),
            ("AcTg", "aCtG"),
            ("actg", "ACTG")
        ]
        
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

    def test_matrix_values(self):
        """Test specific values in the alignment matrix"""
        score, pointer = self.algorithm("AC", "AT")
        
        expected_score = numpy.array([
            [ 0., -5., -6.],
            [-5.,  2., -3.],
            [-6., -3.,  1.]
        ])
        numpy.testing.assert_array_almost_equal(score, expected_score)
        
        expected_pointer = numpy.array([
            [4., 4., 4.],
            [3., 2., 4.],
            [3., 3., 2.]
        ])
        numpy.testing.assert_array_equal(pointer, expected_pointer)

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTG", "ACT"),                  # Longer query
            ("ACT", "ACTG"),                  # Longer subject
            ("ACGT", "AGT"),                  # Internal gap
            ("HOLYROMANEMPIRE", "HOLYPIRE")   # Longer subject
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                score, pointer = self.algorithm(query, subject)
                self.assertEqual(score.shape, (len(query) + 1, len(subject) + 1))
                self.assertTrue(numpy.all(pointer >= 0))  # Valid pointer values

    def test_gap_penalties(self):
        """Test behavior with gaps"""
        test_cases = [
            ("AT", "A-T"),   # Single gap
            ("ACTG", "AG"),  # Multiple gaps
            ("A", "A--A")    # Consecutive gaps
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                _, pointer = self.algorithm(query, subject)
                self.assertTrue(numpy.any(pointer == 3) or numpy.any(pointer == 4))  # Should have gaps

    def test_scoring_parameters(self):
        """Test behavior with different scoring parameters"""
        custom_algorithm = WatermanSmithBeyer(
            match_score=1,
            mismatch_penalty=2,
            new_gap_penalty=3,
            continue_gap_penalty=1
        )
        
        score, _ = custom_algorithm("AC", "AT")
        expected_score = numpy.array([
            [ 0., -4., -5.],
            [-4.,  1., -3.],
            [-5., -3., -1.]
        ])
        numpy.testing.assert_array_almost_equal(score, expected_score)

if __name__ == '__main__':
    unittest.main()
