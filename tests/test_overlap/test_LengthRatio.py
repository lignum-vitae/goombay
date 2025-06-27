import unittest
import numpy
from goombay import LengthRatio


class TestLengthRatio(unittest.TestCase):
    def setUp(self):
        self.algorithm = LengthRatio()

    def test_identical_sequences(self):
        """Test alignment and scoring for identical sequences"""
        seq = "ABC"
        self.assertEqual(self.algorithm.align(seq, seq), f"{seq}\n{seq}")
        self.assertEqual(self.algorithm.similarity(seq, seq), 1.0)
        self.assertEqual(self.algorithm.distance(seq, seq), 0.0)
        self.assertEqual(self.algorithm.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.algorithm.normalized_distance(seq, seq), 0.0)

    def test_completely_different_lengths(self):
        """Test ratios of same character but different lengths"""
        self.assertEqual(self.algorithm.similarity("A", "AAA"), 1 / 3)
        self.assertEqual(self.algorithm.similarity("AAA", "A"), 1 / 3)
        self.assertEqual(self.algorithm.distance("AAA", "A"), 1 - 1 / 3)
        self.assertEqual(self.algorithm.distance("A", "AAA"), 1 - 1 / 3)

    def test_completely_different_lengths_diff(self):
        """Test ratios of same character but different lengths"""
        self.assertEqual(self.algorithm.similarity("T", "AAA"), 1 / 3)
        self.assertEqual(self.algorithm.similarity("AAA", "T"), 1 / 3)
        self.assertEqual(self.algorithm.distance("AAA", "T"), 1 - 1 / 3)
        self.assertEqual(self.algorithm.distance("T", "AAA"), 1 - 1 / 3)

    def test_empty_sequences(self):
        test_cases = [
            ("", "ACTG", 0, 1),
            ("ACTG", "", 0, 1),
            ("", "", 1, 0),
        ]
        for query, subject, sim, dist in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(self.algorithm.distance(query, subject), dist)
                self.assertEqual(
                    self.algorithm.normalized_similarity(query, subject), sim
                )
                self.assertEqual(
                    self.algorithm.normalized_distance(query, subject), dist
                )

    def test_case_sensitivity(self):
        """Checks that algorithm is case insensitive"""
        matrix = self.algorithm.matrix("a", "A")
        self.assertTrue(numpy.array_equal(matrix, numpy.ones((1, 1))))

    def test_matrix_content(self):
        """Matrix should place 1 where characters match"""
        matrix = self.algorithm.matrix("ABC", "AEC")
        expected = numpy.array(
            [
                [1, 0, 0],
                [0, 0, 0],
                [0, 0, 1],
            ]
        )
        numpy.testing.assert_array_equal(matrix, expected)

    def test_align_output(self):
        """Should return query simple alignment of input strings"""
        self.assertEqual(self.algorithm.align("FOO", "BAR"), "FOO\nBAR")

    def test_symmetry(self):
        """Similarity and distance should be symmetric"""
        pairs = [("abc", "abcdef"), ("x", "xyz"), ("short", "longerstring")]
        for query, subject in pairs:
            with self.subTest(query=query, subject=subject):
                self.assertAlmostEqual(
                    self.algorithm.similarity(query, subject),
                    self.algorithm.similarity(subject, query),
                )
                self.assertAlmostEqual(
                    self.algorithm.distance(query, subject),
                    self.algorithm.distance(subject, query),
                )


if __name__ == "__main__":
    unittest.main()
