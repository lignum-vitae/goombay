import unittest
from goombay.algorithms.editdistance import Hamming

class TestHamming(unittest.TestCase):
    """Test suite for Hamming distance algorithm"""
    
    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = Hamming()
        
    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        # Test basic metrics
        self.assertEqual(self.algorithm.similarity("ACTG", "ACTG"), 4)
        self.assertEqual(self.algorithm.distance("ACTG", "ACTG"), 0)
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)
        
    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        # Test basic metrics
        self.assertEqual(self.algorithm.similarity("AAAA", "TTTT"), 0)
        self.assertEqual(self.algorithm.distance("AAAA", "TTTT"), 4)
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.normalized_distance("AAAA", "TTTT"), 1.0)
        
    def test_binary_array(self):
        """Test binary array outputs"""
        test_cases = [
            ("ACTG", "ACTT", [0,0,0,1], [1,1,1,0]),
            ("AAAA", "TTTT", [1,1,1,1], [0,0,0,0]),
            ("ACTG", "ACTG", [0,0,0,0], [1,1,1,1]),
            ("GGGG", "ACTG", [1,1,1,0], [0,0,0,1])
        ]
        
        for query, subject, expected_dist, expected_sim in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.binary_distance_array(query, subject), expected_dist)
                self.assertEqual(self.algorithm.binary_similarity_array(query, subject), expected_sim)
        
    def test_different_lengths(self):
        """Test behavior with different length sequences"""
        test_cases = [
            ("ACTG", "ACT"),
            ("ACT", "ACTG"),
            ("A", "ACTG")
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                with self.assertRaises(ValueError):
                    self.algorithm.distance(query, subject)
                with self.assertRaises(ValueError):
                    self.algorithm.similarity(query, subject)
            
    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "ACTG"),
            ("ACTG", ""),
            ("", "")
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                if query == subject == "":
                    self.assertEqual(self.algorithm.distance(query, subject), 0)
                    self.assertEqual(self.algorithm.similarity(query, subject), 1)
                else:
                    with self.assertRaises(ValueError):
                        self.algorithm.distance(query, subject)
                    with self.assertRaises(ValueError):
                        self.algorithm.similarity(query, subject)
        
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
                    self.algorithm.distance(query, subject),
                    self.algorithm.distance(query.upper(), subject.upper())
                )
                self.assertEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(query.upper(), subject.upper())
                )
        
    def test_integer_input(self):
        """Test behavior with integer inputs"""
        test_cases = [
            (15, 5, 2, 2),    # 1111 vs 0101
            (15, 13, 1, 3),   # 1111 vs 1101
            (15, 15, 0, 4),   # 1111 vs 1111
            (15, 1, 3, 1),    # 1111 vs 0001
            (15, 0, 4, 0)     # 1111 vs 0000
        ]
        
        for num1, num2, exp_dist, exp_sim in test_cases:
            with self.subTest(num1=num1, num2=num2):
                self.assertEqual(self.algorithm.distance(num1, num2), exp_dist)
                self.assertEqual(self.algorithm.similarity(num1, num2), exp_sim)
        
    def test_mixed_input(self):
        """Test behavior with mixed string/int inputs"""
        test_cases = [
            ("ACTG", 5),
            (15, "ACTG")
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                with self.assertRaises(TypeError):
                    self.algorithm.distance(query, subject)
                with self.assertRaises(TypeError):
                    self.algorithm.similarity(query, subject)
        
    def test_call_method(self):
        """Test __call__ method output"""
        test_cases = [
            ("ACTG", "ACTT", 1, [0,0,0,1]),
            (15, 5, 2, [1,0,1,0])
        ]
        
        for query, subject, exp_dist, exp_array in test_cases:
            with self.subTest(query=query, subject=subject):
                dist, dist_array = self.algorithm(query, subject)
                self.assertEqual(dist, exp_dist)
                self.assertEqual(dist_array, exp_array)

    def test_normalization(self):
        """Test normalization behavior"""
        test_cases = [
            ("ACTG", "ACTG", 1.00, 0.00),
            ("ACTG", "AATG", 0.75, 0.25),
            ("ACTG", "AAAG", 0.50, 0.50),
            ("ACTG", "AAAA", 0.25, 0.75), 
            ("ACTG", "VVVV", 0.00, 1.00)
        ]
        
        for query, subject, exp_sim, exp_dist in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.normalized_similarity(query, subject), exp_sim)
                self.assertEqual(self.algorithm.normalized_distance(query, subject), exp_dist)

    def test_alignments(self):
        """Test alignments for Hamming"""
        test_cases = [
            ("ACTG", "ACTG", "ACTG\nACTG"),
            ("A", "A", "A\nA"),
            ("T", "A", "T\nA"),
            ("A", "T", "A\nT"),
            ("", "", "\n"),
            ("TTTT", "AAAA", "TTTT\nAAAA"),
        ]

        for query, subject, alignment in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), alignment)

if __name__ == '__main__':
    unittest.main()
