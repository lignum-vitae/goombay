import unittest
from goombay.algorithms.editdistance import Hirschberg, needleman_wunsch

class TestHirschberg(unittest.TestCase):
    """Test suite for Hirschberg algorithm"""
    
    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = Hirschberg()

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        query = "ACTG"
        subject = "ACTG"
        # Test basic metrics
        self.assertEqual(self.algorithm.distance(query, subject), 0.0)
        self.assertEqual(self.algorithm.similarity(query, subject), 4 * self.algorithm.match_score)
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity(query, subject), 1.0)
        self.assertEqual(self.algorithm.normalized_distance(query, subject), 0.0)
        
        # Test alignment
        self.assertEqual(self.algorithm.align(query, subject), f"{query}\n{subject}")

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        query = "ACTG"
        subject = "FHYU"
        # Test basic metrics
        self.assertEqual(self.algorithm.distance(query, subject), 4*self.algorithm.mismatch_penalty)
        self.assertEqual(self.algorithm.similarity(query, subject), 0.0)
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.algorithm.normalized_distance(query, subject), 1.0)

        #Test alignment
        self.assertEqual(self.algorithm.align(query, subject), f"{query}\n{subject}")

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "ABC"),     # Empty query
            ("ABC", ""),     # Empty subject
            ("", "")         # Both empty
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                h_align = self.algorithm.align(query, subject)
                nw_align = needleman_wunsch.align(query, subject)
                self.assertEqual(h_align, nw_align)
                
                if query == subject == "":
                    self.assertEqual(self.algorithm.similarity(query, subject), 1)
                    self.assertEqual(self.algorithm.distance(query, subject), 0.0)
                else:
                    self.assertEqual(self.algorithm.similarity(query, subject), 0.0)
                    self.assertEqual(self.algorithm.distance(query, subject), self.algorithm.gap_penalty * len(query or subject))

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.align("A", "A"), "A\nA")
        self.assertEqual(self.algorithm.similarity("A", "A"), 1 * self.algorithm.match_score)
        self.assertEqual(self.algorithm.distance("A", "A"), 0.0)
        
        # Test mismatch
        self.assertEqual(self.algorithm.align("A", "T"), "A\nT")
        self.assertEqual(self.algorithm.similarity("A", "T"), 0.0)
        self.assertEqual(self.algorithm.distance("A", "T"), 1 * self.algorithm.mismatch_penalty)

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

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTG", "ACT"),      # Query longer
            ("ACT", "ACTG"),      # Subject longer
            ("AAAAA", "A"),       # Much longer query
            ("A", "AAAAA")        # Much longer subject
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                h_align = self.algorithm.align(query, subject)
                nw_align = needleman_wunsch.align(query, subject)
                self.assertEqual(h_align, nw_align)

    def test_compare_with_needleman(self):
        """Compare alignments with Needleman-Wunsch"""
        test_cases = [
            ("WATERFALL", "WATERFOWL"),
            ("HOLYWATERISABLESSING", "HWATISBLESSING"),
            ("THISISATEST", "THISISNOTASPLENDIDTEST"),
            ("AAAAGGGG", "AAAACCCCGGGG")
        ]
        
        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                h_align = self.algorithm.align(query, subject)
                nw_align = needleman_wunsch.align(query, subject)
                self.assertEqual(h_align, nw_align)

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
