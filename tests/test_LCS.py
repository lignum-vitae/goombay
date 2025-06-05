import unittest
import numpy
from goombay import LongestCommonSubsequence

class TestLCS(unittest.TestCase):
    """Test suite for Longest Common Subsequence algorithm"""
    
    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = LongestCommonSubsequence()
        
    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        # Test alignment
        self.assertEqual(self.algorithm.align("ACTG", "ACTG"), "ACTG")
        
        # Test matrix values
        matrix = self.algorithm("ACTG", "ACTG")
        self.assertEqual(matrix[-1,-1], 4)  # Should find all 4 matches
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)
        
    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        # Test alignment
        self.assertEqual(self.algorithm.align("AAAA", "TTTT"), "There is no common subsequence!")
        
        # Test matrix values
        matrix = self.algorithm("AAAA", "TTTT")
        self.assertEqual(matrix[-1,-1], 0)  # Should find no matches
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.normalized_distance("AAAA", "TTTT"), 1.0)
        
    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "ACTG", 4),    # Empty query
            ("ACTG", "", 4),    # Empty subject
            ("", "", 0)         # Both empty
        ]
        
        for query, subject, dist in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.align(query, subject),
                    "There is no common subsequence!"
                )
                self.assertEqual(self.algorithm.similarity(query, subject), 0)
                self.assertEqual(self.algorithm.distance(query, subject), dist)
                matrix = self.algorithm(query, subject)
                self.assertEqual(matrix[-1,-1], 0)
        
    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.align("A", "A"), "A")
        matrix = self.algorithm("A", "A")
        self.assertEqual(matrix[-1,-1], 1)
        
        # Test mismatch
        self.assertEqual(
            self.algorithm.align("A", "T"),
            "There is no common subsequence!"
        )
        matrix = self.algorithm("A", "T")
        self.assertEqual(matrix[-1,-1], 0)
        
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
                matrix1 = self.algorithm(query, subject)
                matrix2 = self.algorithm(query.upper(), subject.upper())
                numpy.testing.assert_array_equal(matrix1, matrix2)
        
    def test_matrix_shape(self):
        """Test matrix dimensions"""
        query = "ACTG"
        subject = "ACT"
        matrix = self.algorithm(query, subject)
        self.assertEqual(
            matrix.shape,
            (len(query) + 1, len(subject) + 1)
        )
        
    def test_matrix_values(self):
        """Test specific values in the alignment matrix"""
        matrix = self.algorithm("AC", "AT")
        expected = numpy.array([
            [0, 0, 0],
            [0, 1, 1],
            [0, 1, 1]
        ])
        numpy.testing.assert_array_equal(matrix, expected)
        
    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTG", "ACT", "ACT"),          # Longer query
            ("ACT", "ACTG", "ACT"),          # Longer subject
            ("ABCDE", "ACE", "ACE"),         # Internal gaps
            ("HUMAN", "CHIMPANZEE", "HMAN")  # Complex case
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)
                
    def test_repeated_characters(self):
        """Test behavior with repeated characters"""
        test_cases = [
            ("AAAAAA", "AAA", "AAA"),
            ("ABABAB", "ABAB", "ABAB"),
            ("AAAAAA", "TTTTTT", "There is no common subsequence!")
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

if __name__ == '__main__':
    unittest.main()
