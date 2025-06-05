import unittest
import numpy
from goombay import ShortestCommonSupersequence

class TestSCS(unittest.TestCase):
    """Test suite for Shortest Common Supersequence algorithm"""
    
    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = ShortestCommonSupersequence()

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "", "", 0, 0),                # Both empty
            ("", "ACTG", "ACTG", 4, 0),        # Empty query
            ("ACTG", "", "ACTG", 4, 0)        # Empty subject
        ]
        
        for query, subject, expected, dist, sim in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)
                self.assertEqual(self.algorithm.distance(query, subject), dist)
                self.assertEqual(self.algorithm.similarity(query, subject), sim)

    def test_single_character(self):
        """Test behavior with single character sequences"""
        test_cases = [
            ("A", "A", "A"),             # Same character
            ("A", "T", ["AT", "TA"]),    # Different characters
            ("T", "A", ["AT", "TA"])     # Different characters, reversed
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                result = self.algorithm.align(query, subject)
                self.assertIn(result, expected)

    def test_dna_sequences(self):
        """Test behavior with DNA sequences"""
        test_cases = [
            ("ACTG", "ACTG", "ACTG"),            # Identical
            ("AAAA", "TTTT", "AAAATTTT"),        # No common subsequence
            ("AACTG", "ACTGG", "AACTGG"),        # Nested overlap
            ("AGCT", "TACG", "TAGCTG"),          # Partial overlap
            ("AGCT", "GC", "AGCT"),              # Complete overlap
            ("CCTT", "GGGA", "CCTTGGGA")         # No overlap
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                result = self.algorithm.align(query, subject)
                self.assertEqual(len(result), len(expected))
                self.assertEqual(result, expected)

    def test_case_sensitivity(self):
        """Test that matching is case-insensitive"""
        test_cases = [
            ("ACTG", "actg", "ACTG"),
            ("AcTg", "aCtG", "ACTG"),
            ("actg", "ACTG", "ACTG")
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

    def test_matrix(self):
        """Test matrix calculation"""
        test_cases = [
            # Empty sequences
            ("", "", numpy.zeros((1, 1))),
            # Single character
            ("A", "A", numpy.array([
                                   [0, 1], 
                                   [1, 0]
                                   ])),
            # Simple sequence
            ("AT", "TA", numpy.array([
                                     [0, 1, 2], 
                                     [1, 2, 1], 
                                     [2, 1, 2]
                                     ]))
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                result = self.algorithm.matrix(query, subject)
                numpy.testing.assert_array_equal(result, expected)

    def test_similarity(self):
        """Test similarity calculation"""
        test_cases = [
            ("ACTG", "ACTG", 4.0),        # Identical - all positions match
            ("AAAA", "TTTT", 0.0),        # No common subsequence
            ("AACTG", "ACTGG", 4.0),      # ACTG is common subsequence
            ("", "", 0.0),                # Empty sequences
            ("ACTG", "", 0.0),            # One empty sequence
            ("A", "T", 0.0),              # No matching positions in AT/TA
            ("GGCC", "CCGG", 2.0)         # GC is common subsequence
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.similarity(query, subject), expected)

    def test_distance(self):
        """Test distance calculation"""
        test_cases = [
            ("ACTG", "ACTG", 0),         # Identical
            ("AAAA", "TTTT", 8),         # No common subsequence
            ("AACTG", "ACTGG", 2),       # Nested overlap
            ("", "", 0),                 # Empty sequences
            ("ACTG", "", 4),             # One empty sequence
            ("A", "T", 2),               # Single different characters
            ("AATCG", "ATCGG", 2)        # Example from prompt
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.distance(query, subject), expected)

    def test_normalized_similarity(self):
        """Test normalized similarity calculation"""
        test_cases = [
            ("ACTG", "ACTG", 1.0),       # Identical
            ("AAAA", "TTTT", 0.0),       # No common subsequence
            ("AACTG", "ACTGG", 0.667),   # Nested overlap
            ("", "", 1.0),               # Empty sequences
            ("ACTG", "", 0.0),           # One empty sequence
            ("A", "T", 0.0)              # Single different characters
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.normalized_similarity(query, subject),
                    expected,
                    places=3
                )

    def test_normalized_distance(self):
        """Test normalized distance calculation"""
        test_cases = [
            ("ACTG", "ACTG", 0.0),       # Identical
            ("AAAA", "TTTT", 1.0),       # No common subsequence
            ("AACTG", "ACTGG", 0.333),   # Nested overlap
            ("", "", 0.0),               # Empty sequences
            ("ACTG", "", 1.0),           # One empty sequence
            ("A", "T", 1.0)              # Single different characters
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.normalized_distance(query, subject),
                    expected,
                    places=3
                )

if __name__ == '__main__':
    unittest.main()
