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
        self.assertEqual(algo.is_similar("ABC", "ABD"), 0)

    def test_case_insensitivity(self):
        self.assertEqual(self.algorithm.similarity("abc", "ABC"), 0.0)

    def test_empty_sequences(self):
        self.assertEqual(self.algorithm.similarity("", ""), 0.0)
        self.assertEqual(self.algorithm.distance("", ""), 1.0)
        self.assertEqual(self.algorithm.is_similar("", ""), True)

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
