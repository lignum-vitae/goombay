import unittest
import numpy
from goombay import Postfix


class TestPostfix(unittest.TestCase):
    def setUp(self):
        self.algorithm = Postfix()

    def test_identical_sequences(self):
        """Should return full match for identical sequences"""
        seq = "ACTG"
        self.assertEqual(self.algorithm.similarity(seq, seq), len(seq))
        self.assertEqual(self.algorithm.distance(seq, seq), 0)
        self.assertEqual(self.algorithm.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.algorithm.normalized_distance(seq, seq), 0.0)
        self.assertEqual(self.algorithm.align(seq, seq), seq)

    def test_partial_suffix_match(self):
        """Only part of the suffix matches"""
        query = "TTGAC"
        subject = "GGGAC"
        self.assertEqual(self.algorithm.similarity(query, subject), 3)
        self.assertEqual(self.algorithm.distance(query, subject), 2)
        self.assertAlmostEqual(
            self.algorithm.normalized_similarity(query, subject), 0.6
        )
        self.assertAlmostEqual(self.algorithm.normalized_distance(query, subject), 0.4)
        self.assertEqual(self.algorithm.align(query, subject), "GAC")

    def test_no_match(self):
        """No matching suffix at all"""
        query = "ACTG"
        subject = "TGCA"
        self.assertEqual(self.algorithm.similarity(query, subject), 0)
        self.assertEqual(self.algorithm.distance(query, subject), 4)
        self.assertEqual(self.algorithm.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.algorithm.normalized_distance(query, subject), 1.0)
        self.assertEqual(self.algorithm.align(query, subject), "")

    def test_empty_sequences(self):
        """Edge cases with empty input"""
        self.assertEqual(self.algorithm.similarity("", ""), 0)
        self.assertEqual(self.algorithm.distance("", ""), 0)
        with self.assertRaises(ValueError):
            self.algorithm.normalized_similarity("", "")
        with self.assertRaises(ValueError):
            self.algorithm.normalized_distance("", "")
        self.assertEqual(self.algorithm.align("", ""), "")

    def test_single_character(self):
        """Test with one-character sequences"""
        self.assertEqual(self.algorithm.similarity("A", "A"), 1)
        self.assertEqual(self.algorithm.similarity("A", "T"), 0)
        self.assertEqual(self.algorithm.align("A", "A"), "A")
        self.assertEqual(self.algorithm.align("A", "T"), "")

    def test_case_insensitivity(self):
        """Suffix matching should be case-insensitive"""
        self.assertEqual(self.algorithm.similarity("actg", "TACTG"), 4)
        self.assertEqual(self.algorithm.align("aCtG", "TACTG"), "ACTG")

    def test_matrix_output(self):
        """Test matrix for suffix match positions (on diagonals from end)"""
        query = "TTACG"
        subject = "GGACG"
        result = self.algorithm.matrix(query, subject)
        expected = numpy.zeros((5, 5))
        expected[0, 0] = 1
        expected[1, 1] = 1
        expected[2, 2] = 1
        numpy.testing.assert_array_equal(result, expected)

    def test_asymmetric_lengths(self):
        """Suffix match should work even with unequal lengths"""
        query = "CG"
        subject = "ATACG"
        self.assertEqual(self.algorithm.similarity(query, subject), 2)
        self.assertEqual(self.algorithm.align(query, subject), "CG")
        result = self.algorithm.matrix(query, subject)
        expected = numpy.zeros((2, 5))
        expected[0, 0] = 1
        expected[1, 1] = 1
        numpy.testing.assert_array_equal(result, expected)

    def test_suffix_vs_full_match(self):
        """Tests characters are not matched from beginning or middle"""
        query = "TTT"
        subject = "TTTT"
        self.assertEqual(self.algorithm.similarity(query, subject), 3)
        self.assertEqual(self.algorithm.align(query, subject), "TTT")

        query = "TTT"
        subject = "ATTT"
        self.assertEqual(self.algorithm.similarity(query, subject), 3)
        self.assertEqual(self.algorithm.align(query, subject), "TTT")


if __name__ == "__main__":
    unittest.main()
