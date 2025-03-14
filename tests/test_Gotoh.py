import unittest
import numpy
from goombay import Gotoh

class TestGotoh(unittest.TestCase):
    """Test suite for Gotoh global alignment algorithm"""
    
    def setUp(self):
        """Initialize algorithm for tests"""
        self.algorithm = Gotoh()

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
        expected_score = -(4 * self.algorithm.mismatch_penalty)  # All positions are mismatches
        self.assertEqual(self.algorithm.similarity("AAAA", "TTTT"), expected_score)
        
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
                    self.assertEqual(self.algorithm.similarity(query, subject), self.algorithm.match_score)
                    self.assertEqual(self.algorithm.distance(query, subject), 0.0)
                else:
                    non_empty = query or subject
                    expected_gaps = "-" * len(non_empty)
                    self.assertEqual(
                        self.algorithm.align(query, subject),
                        f"{expected_gaps if not query else query}\n{expected_gaps if not subject else subject}"
                    )
                    expected_score = -(self.algorithm.new_gap_penalty + 
                                    len(non_empty) * self.algorithm.continue_gap_penalty)
                    self.assertEqual(self.algorithm.similarity(query, subject), expected_score)

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.align("A", "A"), "A\nA")
        self.assertEqual(self.algorithm.similarity("A", "A"), self.algorithm.match_score)
        self.assertEqual(self.algorithm.distance("A", "A"), 0.0)
        
        # Test mismatch
        self.assertEqual(self.algorithm.align("A", "T"), "A\nT")
        expected_score = -self.algorithm.mismatch_penalty
        self.assertEqual(self.algorithm.similarity("A", "T"), expected_score)

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
        D, P, Q, pointer = self.algorithm(query, subject)
        expected_shape = (len(query) + 1, len(subject) + 1)
        
        self.assertEqual(D.shape, expected_shape)
        self.assertEqual(P.shape, expected_shape)
        self.assertEqual(Q.shape, expected_shape)
        self.assertEqual(pointer.shape, expected_shape)

    def test_gap_extension(self):
        """Test gap extension behavior"""
        query = "ACGGCT"
        # Score with one gap of length 3
        subject = "ACT"
        score1 = self.algorithm.similarity(query, subject)
        self.assertEqual(self.algorithm.align(query, subject), "ACGGCT\nAC---T")
        
        # Score with two gaps of length 2 and 1 respectively
        subject2 = "AGT"
        score2 = self.algorithm.similarity(query, subject2)
        self.assertEqual(self.algorithm.align(query, subject2), "ACGGCT\nA--G-T")
        
        # One long gap should cost less than two short gaps
        self.assertGreater(score1, score2)

    def test_scoring_parameters(self):
        """Test behavior with different scoring parameters"""
        custom_algorithm = Gotoh(
            match_score=1,
            mismatch_penalty=2,
            new_gap_penalty=3,
            continue_gap_penalty=1
        )
        
        # Test alignment
        self.assertEqual(custom_algorithm.align("ACGT", "AGT"), "ACGT\nA-GT")
        
        # Test that gap extension is cheaper than gap opening
        query = "ACGT"
        subject = "AT"
        D, _, _, _ = custom_algorithm(query, subject)
        expected_penalty = (custom_algorithm.new_gap_penalty + 2*custom_algorithm.continue_gap_penalty)
        expected_matches = 2 * custom_algorithm.match_score
        self.assertEqual(D[-1, -1], expected_matches - expected_penalty)

    def test_matrix_values(self):
        """Test specific matrix values"""
        D, P, Q, pointer = self.algorithm("AC", "AT")
        
        # First row and column should reflect gap penalties
        self.assertTrue(numpy.all(D[0, 1:] <= 0))  # First row
        self.assertTrue(numpy.all(D[1:, 0] <= 0))  # First column
        
        # Gap matrices initialization
        self.assertTrue(numpy.all(P[:, 0] == 0))  # First column of P should be 0
        self.assertTrue(numpy.all(Q[0, :] == 0))  # First row of Q should be 0
        self.assertTrue(numpy.all(P[:, 1:] <= 0))  # Rest of P should be non-positive
        self.assertTrue(numpy.all(Q[1:, :] <= 0))  # Rest of Q should be non-positive
        
        # Test pointer values are valid
        self.assertTrue(numpy.all(pointer >= 0))
        self.assertTrue(numpy.all(pointer <= 9))  # Max valid pointer value

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
