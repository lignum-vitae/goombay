import unittest
import numpy
from goombay import MLIPNS


class TestMLIPNS(unittest.TestCase):
    def setUp(self):
        self.algorithm = MLIPNS()

    def test_identical_sequences(self):
        seq = "PRODUCT"
        self.assertEqual(self.algorithm.align(seq, seq), f"{seq}\n{seq}")
        self.assertEqual(self.algorithm.similarity(seq, seq), 0.0)
        self.assertEqual(self.algorithm.distance(seq, seq), 1.0)
        self.assertEqual(self.algorithm.is_similar(seq, seq), True)

    def test_mismatch_within_limit(self):
        self.assertEqual(self.algorithm.is_similar("ABC", "ABD"), True)
        matrix = self.algorithm("ABC", "ABD")
        self.assertTrue(numpy.sum(matrix) > 0)

    def test_mismatch_exceeds_limit(self):
        algo = MLIPNS(max_mismatch=0)
        matrix = algo("ABC", "ABD")
        self.assertTrue(numpy.all(matrix == 0))
        self.assertEqual(algo.is_similar("ABC", "ABD"), False)

    def test_case_insensitivity(self):
        self.assertEqual(self.algorithm.similarity("abc", "ABC"), 0.0)

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

    def test_symmetry(self):
        pairs = [("abc", "cab"), ("Tomato", "Tamato"), ("ABCD", "AXYZ")]
        for query, subject in pairs:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(subject, query),
                )

    def test_number_deletions(self):
        """Check is_similar behavior under different threshold settings in MLIPNS"""
        # 2 or fewer deletions
        query, subject = "ABCD", "ABYZ"
        self.assertEqual(self.algorithm.is_similar(query, subject), True)

        # more than 2 deletions
        query, subject = "ABCD", "AXYZ"
        self.assertEqual(self.algorithm.is_similar(query, subject), False)


if __name__ == "__main__":
    unittest.main()
