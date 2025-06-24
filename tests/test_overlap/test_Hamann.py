import unittest
import numpy
from goombay import Hamann


class TestHamann(unittest.TestCase):
    def setUp(self):
        self.algorithm = Hamann()

    def test_identical_strings(self):
        """Test scoring of identical sequence"""
        seq = "ACGT"
        self.assertAlmostEqual(self.algorithm.similarity(seq, seq), 1.0)
        self.assertAlmostEqual(self.algorithm.distance(seq, seq), 0.0)
        self.assertAlmostEqual(self.algorithm.normalized_similarity(seq, seq), 1.0)
        self.assertAlmostEqual(self.algorithm.normalized_distance(seq, seq), 0.0)

    def test_half_match(self):
        """Two matching, two mismatching characters"""
        query = "AAGG"
        subject = "AATT"
        self.assertAlmostEqual(self.algorithm.similarity(query, subject), 0.0)
        self.assertAlmostEqual(
            self.algorithm.normalized_similarity(query, subject), 0.5
        )
        self.assertAlmostEqual(self.algorithm.distance(query, subject), 0.5)
        self.assertAlmostEqual(self.algorithm.normalized_distance(query, subject), 0.5)

    def test_completely_different(self):
        """Test no characters matching"""
        query = "AAAA"
        subject = "TTTT"
        self.assertAlmostEqual(self.algorithm.similarity(query, subject), -1.0)
        self.assertAlmostEqual(
            self.algorithm.normalized_similarity(query, subject), 0.0
        )
        self.assertAlmostEqual(self.algorithm.distance(query, subject), 1.0)
        self.assertAlmostEqual(self.algorithm.normalized_distance(query, subject), 1.0)

    def test_binary(self):
        """Check binary co-occurrence matrix"""
        query = "1010"
        subject = "1100"
        expected = numpy.array(
            [
                [1, 1],
                [1, 1],
            ]
        )
        result = self.algorithm.matrix(query, subject, binary=True)
        numpy.testing.assert_array_equal(result, expected)

    def test_character_matrix(self):
        """Check character match/mismatch matrix"""
        query = "ACGT"
        subject = "AGGT"
        result = self.algorithm.matrix(query, subject, binary=False)
        expected = numpy.array([3, 1])
        self.assertEqual(result.shape, (2,))
        numpy.testing.assert_array_equal(result, expected)

    def test_input_validation(self):
        """Should raise errors on invalid input types or unequal lengths"""
        with self.assertRaises(TypeError):
            self.algorithm.similarity(["A", "C"], "AC")
        with self.assertRaises(ValueError):
            self.algorithm.similarity("AC", "ACGT")

    def test_symmetry(self):
        """Hamann similarity should be symmetric"""
        a = "ACGT"
        b = "AGGT"
        sim1 = self.algorithm.similarity(a, b)
        sim2 = self.algorithm.similarity(b, a)
        self.assertAlmostEqual(sim1, sim2)


if __name__ == "__main__":
    unittest.main()
