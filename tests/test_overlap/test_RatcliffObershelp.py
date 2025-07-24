import unittest
from goombay import RatcliffObershelp


class TestRatcliffObershelp(unittest.TestCase):
    """Test suite for the Ratcliff‑Obershelp similarity algorithm"""

    def setUp(self):
        self.algorithm = RatcliffObershelp()

    def test_identical_sequences(self):
        """Identical inputs should match completely"""
        self.assertEqual(self.algorithm.align("ACTG", "ACTG"), ["ACTG"])
        self.assertEqual(self.algorithm.similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.distance("ACTG", "ACTG"), 0.0)
        self.assertEqual(self.algorithm.normalized_similarity("ACTG", "ACTG"), 1.0)
        self.assertEqual(self.algorithm.normalized_distance("ACTG", "ACTG"), 0.0)

    def test_completely_different(self):
        """No common characters → empty alignment"""
        self.assertEqual(self.algorithm.align("AAAA", "TTTT"), [])
        self.assertEqual(self.algorithm.similarity("AAAA", "TTTT"), 0.0)
        self.assertEqual(self.algorithm.distance("AAAA", "TTTT"), 1.0)

    def test_empty_sequences(self):
        """All combinations of empty inputs"""
        test_cases = [
            ("", "ACTG", 0.0, 1.0),  # Empty query
            ("ACTG", "", 0.0, 1.0),  # Empty subject
            ("", "", 1.0, 0.0),  # Both empty
        ]

        for query, subject, sim, dist in test_cases:
            with self.subTest(query=query, subject=subject):
                self.assertEqual(self.algorithm.align(query, subject), [])
                self.assertEqual(self.algorithm.similarity(query, subject), sim)
                self.assertEqual(self.algorithm.distance(query, subject), dist)

    def test_single_character(self):
        """Match and mismatch for one‑character strings"""
        # Perfect match
        self.assertEqual(self.algorithm.align("A", "A"), ["A"])
        self.assertEqual(self.algorithm.similarity("A", "A"), 1.0)
        # Mismatch
        self.assertEqual(self.algorithm.align("A", "T"), [])
        self.assertEqual(self.algorithm.similarity("A", "T"), 0.0)

    def test_case_sensitivity(self):
        """The algorithm should treat upper‑ and lower‑case the same"""
        test_pairs = [("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")]

        for q, s in test_pairs:
            with self.subTest(query=q, subject=s):
                self.assertEqual(
                    self.algorithm.align(q, s),
                    self.algorithm.align(q.upper(), s.upper()),
                )
                self.assertEqual(
                    self.algorithm.similarity(q, s),
                    self.algorithm.similarity(q.upper(), s.upper()),
                )

    def test_different_lengths(self):
        """Typical L‑shape alignments where lengths differ"""
        cases = [
            ("ACTG", "ACT", ["ACT"]),  # Extra char at end of query
            ("ACT", "ACTG", ["ACT"]),  # Extra char at end of subject
            ("ABCDE", "ACE", ["A", "C", "E"]),  # Internal gaps
            ("HUMAN", "CHIMPANZEE", ["AN", "H", "M"]),  # More complex
        ]
        for q, s, expected in cases:
            with self.subTest(query=q, subject=s):
                self.assertEqual(self.algorithm.align(q, s), expected)

    def test_repeated_characters(self):
        """Ensure repeated runs are handled correctly"""
        cases = [
            ("AAAAAA", "AAA", ["AAA"]),
            ("ABABAB", "ABAB", ["ABAB"]),
            ("AAAAAA", "TTTTTT", []),
        ]
        for q, s, expected in cases:
            with self.subTest(query=q, subject=s):
                self.assertEqual(self.algorithm.align(q, s), expected)

    def test_symmetry(self):
        """similarity(a, b) must equal similarity(b, a)"""
        pairs = [
            ("KITTEN", "SITTING"),
            ("FLIGHT", "NIGHT"),
            ("ABCDE", "EDCBA"),
        ]
        for a, b in pairs:
            with self.subTest(a=a, b=b):
                self.assertEqual(
                    self.algorithm.similarity(a, b),
                    self.algorithm.similarity(b, a),
                )
                self.assertEqual(
                    self.algorithm.distance(a, b),
                    self.algorithm.distance(b, a),
                )

    def test_known_values(self):
        cases = [
            ("WIKIMEDIA", "WIKIMANIA", (2 * (5 + 2)) / (9 + 9)),
            ("GESTALT PATTERN MATCHING", "GESTALT PRACTICE", (2 * 13) / (24 + 16)),
        ]
        for q, s, sim in cases:
            with self.subTest(query=q, subject=s):
                self.assertEqual(self.algorithm.similarity(q, s), sim)


if __name__ == "__main__":
    unittest.main()
