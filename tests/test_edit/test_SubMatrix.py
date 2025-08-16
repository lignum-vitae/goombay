import unittest
from biobase.matrix import Blosum, Pam
from goombay import NeedlemanWunsch, Gotoh, Hirschberg, WatermanSmithBeyer


class TestSubstitutionMatrices(unittest.TestCase):
    """Test suite for substitution matrices that can be used with global alignment algorithms"""

    def setUp(self):
        """Initialize algorithm for tests"""
        self.nw = NeedlemanWunsch()
        self.g = Gotoh()
        self.h = Hirschberg()
        self.nwb62 = NeedlemanWunsch(scoring_matrix=Blosum(62))
        self.hb62 = Hirschberg(scoring_matrix=Blosum(62))
        self.wsbp250 = WatermanSmithBeyer(
            new_gap=5, continued_gap=3, scoring_matrix=Pam(250)
        )
        self.gp250 = Gotoh(new_gap=5, continued_gap=3, scoring_matrix=Pam(250))

    def test_identical_sequences(self):
        """Test behavior with identical sequences"""
        seq = "ARLP"

        self.assertEqual(self.nw.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.nw.normalized_distance(seq, seq), 0.0)

        self.assertEqual(self.g.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.g.normalized_distance(seq, seq), 0.0)

        self.assertEqual(self.h.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.h.normalized_distance(seq, seq), 0.0)

        self.assertEqual(self.nwb62.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.nwb62.normalized_distance(seq, seq), 0.0)

        self.assertEqual(self.hb62.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.hb62.normalized_distance(seq, seq), 0.0)

        self.assertEqual(self.wsbp250.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.wsbp250.normalized_distance(seq, seq), 0.0)

        self.assertEqual(self.gp250.normalized_similarity(seq, seq), 1.0)
        self.assertEqual(self.gp250.normalized_distance(seq, seq), 0.0)

    def test_worst_alignment_score(self):
        """Test behavior with identical sequences"""
        query = "*****"
        subject = "DDDDD"

        # Test normalization
        self.assertEqual(self.nw.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.nw.normalized_distance(query, subject), 1.0)

        self.assertEqual(self.g.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.g.normalized_distance(query, subject), 1.0)

        self.assertEqual(self.h.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.h.normalized_distance(query, subject), 1.0)

        self.assertEqual(self.nwb62.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.nwb62.normalized_distance(query, subject), 1.0)

        self.assertEqual(self.hb62.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.hb62.normalized_distance(query, subject), 1.0)

        self.assertEqual(self.wsbp250.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.wsbp250.normalized_distance(query, subject), 1.0)

        self.assertEqual(self.gp250.normalized_similarity(query, subject), 0.0)
        self.assertEqual(self.gp250.normalized_distance(query, subject), 1.0)

    def test_different_length_sequences(self):
        query = "MKT"
        subject = "MKTTT"

        # With identical prefix but gaps at the end
        self.assertGreaterEqual(self.nw.normalized_similarity(query, subject), 0.0)
        self.assertLessEqual(self.nw.normalized_similarity(query, subject), 1.0)

        self.assertGreaterEqual(self.g.normalized_similarity(query, subject), 0.0)
        self.assertLessEqual(self.g.normalized_similarity(query, subject), 1.0)

        self.assertGreaterEqual(self.h.normalized_similarity(query, subject), 0.0)
        self.assertLessEqual(self.h.normalized_similarity(query, subject), 1.0)

        self.assertGreaterEqual(self.nwb62.normalized_similarity(query, subject), 0.0)
        self.assertLessEqual(self.nwb62.normalized_similarity(query, subject), 1.0)

        self.assertGreaterEqual(self.hb62.normalized_similarity(query, subject), 0.0)
        self.assertLessEqual(self.hb62.normalized_similarity(query, subject), 1.0)

        self.assertGreaterEqual(self.wsbp250.normalized_similarity(query, subject), 0.0)
        self.assertLessEqual(self.wsbp250.normalized_similarity(query, subject), 1.0)

    def test_pam250_negative_scores(self):
        query = "AW"
        subject = "WA"
        sim_score = self.wsbp250.normalized_similarity(query, subject)
        dist_score = self.wsbp250.normalized_distance(query, subject)
        self.assertGreaterEqual(sim_score, 0.0)
        self.assertGreaterEqual(dist_score, 0.0)
        self.assertLessEqual(sim_score, 1.0)
        self.assertLessEqual(dist_score, 1.0)

    def test_empty_sequences(self):
        self.assertEqual(self.nw.normalized_similarity("", ""), 1.0)
        self.assertEqual(self.g.normalized_similarity("", ""), 1.0)
        self.assertEqual(self.h.normalized_similarity("", ""), 1.0)

        self.assertEqual(self.nw.normalized_distance("", ""), 0.0)
        self.assertEqual(self.g.normalized_distance("", ""), 0.0)
        self.assertEqual(self.h.normalized_distance("", ""), 0.0)


if __name__ == "__main__":
    unittest.main()
