import unittest
import numpy
from goombay import Prefix


class TestPrefix(unittest.TestCase):
    def setUp(self):
        self.algorithm = Prefix()

    def test_identical_sequences(self):
        """Should return full match for identical sequences"""
        seq = "ACTG"
        self.assertEqual(self.algorithm.similarity(seq, seq), len(seq))
        self.assertEqual(self.algorithm.distance(seq, seq), 0)
        self.assertEqual(self.algorithm.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.algorithm.normalized_distance(seq, seq), 0.0)
        self.assertEqual(self.algorithm.align(seq, seq), seq)

    def test_partial_prefix_match(self):
        """Only part of the prefix matches"""
        query = "ACTG"
        subject = "ACGG"
        self.assertEqual(self.algorithm.similarity(query, subject), 2)
        self.assertEqual(self.algorithm.distance(query, subject), 2)
        self.assertAlmostEqual(
            self.algorithm.normalized_similarity(query, subject), 0.5
        )
        self.assertAlmostEqual(self.algorithm.normalized_distance(query, subject), 0.5)
        self.assertEqual(self.algorithm.align(query, subject), "AC")

    def test_no_match(self):
        """No matching prefix at all"""
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
        self.assertTrue((self.algorithm.matrix("", "") == numpy.zeros((0, 0))).all())

    def test_single_character(self):
        """Test with one-character sequences"""
        self.assertEqual(self.algorithm.similarity("A", "A"), 1)
        self.assertEqual(self.algorithm.similarity("A", "T"), 0)
        self.assertEqual(self.algorithm.align("A", "A"), "A")
        self.assertEqual(self.algorithm.align("A", "T"), "")

    def test_case_sensitivity(self):
        """Prefix matching should be case-insensitive"""
        self.assertEqual(self.algorithm.similarity("actg", "ACTG"), 4)
        self.assertEqual(self.algorithm.similarity("aCt", "ACT"), 3)
        self.assertEqual(self.algorithm.align("aCtG", "ACTG"), "ACTG")

    def test_matrix_output(self):
        """Test that matrix marks prefix match positions with 1s on the diagonal"""
        query = "ACTG"
        subject = "ACAG"
        result = self.algorithm.matrix(query, subject)
        expected = numpy.zeros((4, 4))
        expected[0, 0] = 1
        expected[1, 1] = 1
        numpy.testing.assert_array_equal(result, expected)

    def test_asymmetric_lengths(self):
        """Prefix should match as long as the shorter one allows"""
        query = "ACT"
        subject = "ACTGGA"
        self.assertEqual(self.algorithm.similarity(query, subject), 3)
        self.assertEqual(self.algorithm.align(query, subject), "ACT")
        result = self.algorithm.matrix(query, subject)
        expected = numpy.zeros((3, 6))
        expected[0, 0] = 1
        expected[1, 1] = 1
        expected[2, 2] = 1
        numpy.testing.assert_array_equal(result, expected)

    def test_prefix_vs_full_match(self):
        """Tests that interior matches are ignored"""
        query = "TTT"
        subject = "ATTT"
        self.assertEqual(self.algorithm.similarity(query, subject), 0)
        self.assertEqual(self.algorithm.align(query, subject), "")


if __name__ == "__main__":
    unittest.main()
