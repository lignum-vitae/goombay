import unittest
import numpy
from goombay import Gotoh_Local

class TestGotohLocal(unittest.TestCase):
    """Test suite for Gotoh Local alignment algorithm"""
    
    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = Gotoh_Local()

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        # Test alignment
        self.assertEqual(self.algorithm.align("ACTG", "ACTG"), "ACTG\nACTG")
        
        # Test scoring
        self.assertEqual(self.algorithm.similarity("ACTG", "ACTG"), 4*self.algorithm.match_score)
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)

    def test_completely_different(self):
        """Test behavior with completely different sequences"""
        # Test alignment
        self.assertEqual(self.algorithm.align("AAAA", "TTTT"), "There is no local alignment!")
        
        # Test scoring
        self.assertEqual(self.algorithm.similarity("AAAA", "TTTT"), 0.0)
        
        # Test normalization
        self.assertEqual(self.algorithm.normalized_similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.normalized_distance("AAAA", "TTTT"), 1.0)

    def test_empty_sequences(self):
        """Test behavior with empty sequences"""
        test_cases = [
            ("", "ACTG", 4, 1, 0, 0),    # Empty query
            ("ACTG", "", 4, 1, 0, 0),    # Empty subject
            ("", "", 0, 0, 1, 1),    # Empty subject
        ]
        
        for query, subject, dist, norm_dist, sim, norm_sim in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), "There is no local alignment!")
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(self.algorithm.normalized_similarity(query, subject), norm_sim)
                self.assertEqual(self.algorithm.distance(query, subject), dist)
                self.assertEqual(self.algorithm.normalized_distance(query, subject), norm_dist)

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.align("A", "A"), "A\nA")
        self.assertEqual(self.algorithm.similarity("A", "A"), 1 * self.algorithm.match_score)
        self.assertEqual(self.algorithm.distance("A", "A"), 0.0)
        
        # Test mismatch
        self.assertEqual(self.algorithm.align("A", "T"), "There is no local alignment!")
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
                self.assertEqual(self.algorithm.align(query, subject), "ACTG\nACTG")
                self.assertEqual(self.algorithm.similarity(query, subject), 4*self.algorithm.match_score)
                self.assertEqual(self.algorithm.distance(query, subject), 0.0)

    def test_matrix_shape(self):
        """Test matrix dimensions"""
        query = "ACGT"
        subject = "AGT"
        D, P, Q = self.algorithm.matrix(query, subject)
        expected_shape = (len(query) + 1, len(subject) + 1)
        self.assertEqual(D.shape, expected_shape)
        self.assertEqual(P.shape, expected_shape)
        self.assertEqual(Q.shape, expected_shape)

    def test_matrix_values(self):
        """Test specific matrix values"""
        query = "AC"
        subject = "AG"
        D, P, Q = self.algorithm.matrix(query, subject)
        
        # Test D matrix (main scoring matrix)
        expected_D = numpy.array([
            [0., 0., 0.],
            [0., 2., 0.],
            [0., 0., 1.]
        ])
        numpy.testing.assert_array_equal(D, expected_D)
        
        # Test P and Q matrices (gap matrices)
        self.assertTrue(numpy.all(P <= 0))  # Gap scores should be non-positive
        self.assertTrue(numpy.all(Q <= 0))  # Gap scores should be non-positive

    def test_different_lengths(self):
        """Test behavior with sequences of different lengths"""
        test_cases = [
            ("ACTG", "ACT", "ACT\nACT"),    # Longer query
            ("ACT", "ACTG", "ACT\nACT"),    # Longer subject
            ("ACGT", "AGT", "GT\nGT")       # Internal gap
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

    def test_scoring_parameters(self):
        """Test behavior with different scoring parameters"""
        custom_algorithm = Gotoh_Local(
            match_score=2,
            mismatch_penalty=1,
            new_gap_penalty=2,
            continue_gap_penalty=1
        )
        
        # Test alignment
        self.assertEqual(custom_algorithm.align("ACGT", "AGT"), "GT\nGT")
        self.assertEqual(custom_algorithm.similarity("ACGT", "AGT"), 4.0)

    def test_local_alignment_behavior(self):
        """Test specific local alignment behaviors"""
        test_cases = [
            ("CGATC", "GTATG", "AT\nAT"),              # Find internal match
            ("AAAGGGCCGGTTT", "AAATTT", "AAA\nAAA"),   # Multiple possible alignments
            ("ACGTACGT", "TACGTAC", "ACGTAC\nACGTAC")  # Longer local alignment
        ]
        
        for query, subject, expected in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), expected)

if __name__ == '__main__':
    unittest.main()
