import unittest
from goombay import Jaro_Winkler

class TestJaroWinkler(unittest.TestCase):
    """Test suite for Jaro-Winkler distance algorithm"""
    
    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = Jaro_Winkler()
        
    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        test_cases = [
            "ACTG",          # DNA sequence
            "1010",          # Binary string
            "Hello",         # Regular string
            "A",            # Single character
            ""              # Empty string
        ]
        
        for sequence in test_cases:
            with self.subTest(sequence=sequence):
                # Test similarity
                self.assertEqual(self.algorithm.similarity(sequence, sequence), 1.0)
                
                # Test distance
                self.assertEqual(self.algorithm.distance(sequence, sequence), 0.0)
                
                # Test normalization
                self.assertEqual(self.algorithm.normalized_similarity(sequence, sequence), 1.0)
                self.assertEqual(self.algorithm.normalized_distance(sequence, sequence), 0.0)

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        test_cases = [
            ("AAAA", "TTTT"),          # DNA sequences
            ("0000", "1111"),          # Binary strings
            ("Hello", "Warvd"),        # Regular strings
            ("ABC", "XYZ")             # No matching characters
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                # Test similarity
                self.assertEqual(self.algorithm.similarity(query, subject), 0.0)
                
                # Test distance
                self.assertEqual(self.algorithm.distance(query, subject), 1.0)
                
                # Test normalization
                self.assertEqual(self.algorithm.normalized_similarity(query, subject), 0.0)
                self.assertEqual(self.algorithm.normalized_distance(query, subject), 1.0)

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "", 0, 1, "\n"),             # Both empty sequences
            ("A", "", 1, 0, "A\n-"),          # One empty subject sequence
            ("", "A", 1, 0, "-\nA"),          # One empty query sequence
            ("", "ACTG", 1, 0, "----\nACTG"), # Longer empty query sequence
            ("ACTG", "", 1, 0, "ACTG\n----")  # Longer empty subject sequence
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
        test_cases = [
            ("ACTG", "actg"),
            ("AcTg", "aCtG"),
            ("actg", "ACTG")
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(query.upper(), subject.upper())
                )
                self.assertEqual(
                    self.algorithm.distance(query, subject),
                    self.algorithm.distance(query.upper(), subject.upper())
                )

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTG", "ACT"),      # Longer query
            ("ACT", "ACTG"),      # Longer subject
            ("A", "ACTG"),        # Much shorter query
            ("ACTG", "A"),        # Much shorter subject
            ("T", "ACGTTTGGGAC"), # Significantly shorter query
            ("ACGTTTGGGAC", "T")  # Significantly shorter subject
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                sim = self.algorithm.similarity(query, subject)
                dist = self.algorithm.distance(query, subject)
                self.assertTrue(0 <= sim <= 1)
                self.assertTrue(0 <= dist <= 1)
                self.assertAlmostEqual(sim + dist, 1.0)

    def test_similarity_scores(self):
        """Test specific similarity scores"""
        test_cases = [
            ("FAREMVIEL", "FARMVILLE", 0.919),
            ("CRATE", "TRACE", 0.733),
            ("DWAYNE", "DUANE", 0.840),
            ("DIC", "DICKSON", 0.867),
            ("DIXON", "DICKSONX", 0.813),
            ("MARTHA", "MART", 0.933),
            ("MART", "MARTHA", 0.933),
            ("MARTHA", "MARHTA", 0.961)
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.similarity(query, subject),
                    expected,
                    places=3
                )
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.normalized_similarity(query, subject),
                    expected,
                    places=3
                )

    def test_prefix_bonus(self):
        """Test prefix bonus behavior"""
        # Identical prefixes should get bonus
        self.assertGreater(
            self.algorithm.similarity("MARTHA", "MARHTA"),
            self.algorithm.similarity("ARDWARK", "AARDVARK")
        )
        
        # Length of common prefix affects bonus
        self.assertGreater(
            self.algorithm.similarity("ABCDEF", "ABCXXX"),
            self.algorithm.similarity("ABCDEF", "AXXXXX")
        )

if __name__ == '__main__':
    unittest.main()
