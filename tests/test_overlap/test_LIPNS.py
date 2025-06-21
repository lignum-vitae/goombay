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
        self.assertEqual(self.algorithm.similarity("", ""), 0.0)
        self.assertEqual(self.algorithm.distance("", ""), 1.0)
        self.assertEqual(self.algorithm.is_similar("", ""), True)

    def test_known_solution(self):
        self.assertAlmostEqual(
            self.algorithm.similarity("Tomato", "Tamato"), 0.16, delta=0.01
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
