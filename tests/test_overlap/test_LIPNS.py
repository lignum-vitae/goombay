import unittest
from goombay import LIPNS


class TestLIPNS(unittest.TestCase):
    def setUp(self):
        self.algorithm = LIPNS()

    def test_identical_sequences(self):
        """Should return perfect match for identical sequences"""
        seq = "ABC"
        self.assertEqual(self.algorithm.align(seq, seq), f"{seq}\n{seq}")
        self.assertEqual(self.algorithm.similarity(seq, seq), 0.0)
        self.assertEqual(self.algorithm.distance(seq, seq), 1.0)
        self.assertEqual(self.algorithm.normalized_similarity(seq, seq), 0.0)
        self.assertEqual(self.algorithm.normalized_distance(seq, seq), 1.0)
        self.assertEqual(self.algorithm.is_similar(seq, seq), True)

    def test_case_insensitivity(self):
        self.assertEqual(self.algorithm.similarity("abc", "ABC"), 0.0)

    def test_completely_different(self):
        """No matching characters"""
        self.assertEqual(self.algorithm.similarity("AAA", "BBB"), 1.0)
        self.assertEqual(self.algorithm.normalized_similarity("AAA", "BBB"), 1.0)
        self.assertEqual(self.algorithm.is_similar("AAA", "BBB"), 0)

    def test_empty_sequences(self):
        test_cases = [
            ("", "ACTG", 1, 0, False),
            ("ACTG", "", 1, 0, False),
            ("", "", 0, 1, True),
        ]
        for query, subject, sim, dist, is_sim in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(self.algorithm.distance(query, subject), dist)
                self.assertEqual(self.algorithm.is_similar(query, subject), is_sim)

    def test_single_character(self):
        """Test behavior with single character sequences"""
        # Test match
        self.assertEqual(self.algorithm.similarity("A", "A"), 0)
        self.assertEqual(self.algorithm.distance("A", "A"), 1)
        self.assertEqual(self.algorithm.is_similar("A", "A"), True)

        # Test mismatch
        self.assertEqual(self.algorithm.similarity("A", "T"), 1)
        self.assertEqual(self.algorithm.distance("A", "T"), 0)
        self.assertEqual(self.algorithm.is_similar("A", "T"), False)

    def test_case_sensitivity(self):
        """Test that matching is case-insensitive"""
        test_cases = [("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")]

        for query, subject in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(query.upper(), subject.upper()),
                )
                self.assertEqual(
                    self.algorithm.distance(query, subject),
                    self.algorithm.distance(query.upper(), subject.upper()),
                )

    def test_known_solution(self):
        self.assertAlmostEqual(
            self.algorithm.similarity("Tomato", "Tamato"), 0.16, delta=0.01
        )

    def test_symmetry(self):
        """Checks that algorithm is symmetric"""
        pairs = [("abc", "cab"), ("Tomato", "Tamato"), ("ABCD", "AXYZ")]
        for query, subject in pairs:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(subject, query),
                )

    def test_threshold_behavior(self):
        """Check is_similar behavior under different threshold settings"""
        query, subject = "ABCD", "AXYZ"

        actual_similarity = self.algorithm.similarity(query, subject)
        # High threshold: should allow the match
        high_threshold = LIPNS(threshold=1.0)
        self.assertEqual(high_threshold.is_similar(query, subject), True)

        # Low threshold: should not allow the match
        low_threshold = LIPNS(threshold=actual_similarity - 0.01)
        self.assertEqual(low_threshold.is_similar(query, subject), False)


if __name__ == "__main__":
    unittest.main()
