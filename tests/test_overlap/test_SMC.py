import unittest
from goombay import SimpleMatchingCoefficient


class TestSimpleMatchingCoefficient(unittest.TestCase):
    def setUp(self):
        self.algorithm = SimpleMatchingCoefficient()

    def test_identical_sequences(self):
        """Identical sequences should return perfect similarity and zero distance"""
        seq = "ACGT"
        self.assertAlmostEqual(self.algorithm.similarity(seq, seq), 1.0)
        self.assertAlmostEqual(self.algorithm.distance(seq, seq), 0.0)
        self.assertAlmostEqual(self.algorithm.normalized_similarity(seq, seq), 1.0)
        self.assertAlmostEqual(self.algorithm.normalized_distance(seq, seq), 0.0)

    def test_half_match(self):
        """50% match should yield 0.5 similarity"""
        query = "AAGG"
        subject = "AATT"
        self.assertAlmostEqual(self.algorithm.similarity(query, subject), 0.5)
        self.assertAlmostEqual(self.algorithm.distance(query, subject), 0.5)

    def test_completely_different(self):
        """No matching characters"""
        query = "AAAA"
        subject = "TTTT"
        self.assertAlmostEqual(self.algorithm.similarity(query, subject), 0.0)
        self.assertAlmostEqual(self.algorithm.distance(query, subject), 1.0)

    def test_empty_sequences(self):
        """Empty sequences should return 0.0 similarity and distance"""
        with self.assertRaises(ValueError):
            self.algorithm.similarity("", "")
        with self.assertRaises(ValueError):
            self.algorithm.distance("", "")
        with self.assertRaises(ValueError):
            self.algorithm.normalized_similarity("", "")
        with self.assertRaises(ValueError):
            self.algorithm.normalized_distance("", "")

    def test_case_insensitivity(self):
        """Matching should be case-insensitive"""
        query = "acgt"
        subject = "ACGT"
        self.assertAlmostEqual(self.algorithm.similarity(query, subject), 1.0)

    def test_symmetry(self):
        """Similarity and distance should be symmetric"""
        a = "ACGT"
        b = "AGGT"
        self.assertAlmostEqual(
            self.algorithm.similarity(a, b), self.algorithm.similarity(b, a)
        )
        self.assertAlmostEqual(
            self.algorithm.distance(a, b), self.algorithm.distance(b, a)
        )

    def test_invalid_length(self):
        """Mismatched sequence lengths should raise IndexError"""
        with self.assertRaises(IndexError):
            self.algorithm.similarity("ACG", "ACGT")

    def test_invalid_type(self):
        """Non-string inputs should raise TypeError inside hamming"""
        with self.assertRaises(TypeError):
            self.algorithm.similarity(["A", "C"], ["A", "C"])


if __name__ == "__main__":
    unittest.main()
