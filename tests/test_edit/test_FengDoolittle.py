import unittest
from goombay import FengDoolittle


class TestFengDoolittle(unittest.TestCase):
    """Test suite for Feng-Doolittle progressive multiple sequence alignment"""

    def setUp(self):
        """Initialize algorithm with default pairwise/clustering (NW + NJ)"""
        self.feng = FengDoolittle()

    def test_identical_sequences(self):
        """Align multiple identical sequences"""
        seqs = ["ACTG", "ACTG", "ACTG"]
        result = self.feng.align(seqs).splitlines()

        print(result)
        self.assertEqual(len(result), 3)
        for line in result:
            self.assertEqual(line, "ACTG")

    def test_progressive_alignment_order(self):
        """Check output for a known 3-sequence progressive alignment"""
        seqs = ["HOUSEOFCARDSFALLDOWN", "HOUSECARDFALLDOWN", "FALLDOWN"]
        result = self.feng.align(seqs).splitlines()

        self.assertEqual(len(result), 3)
        self.assertTrue(all(len(row) == len(result[0]) for row in result))  # aligned
        self.assertIn("FALLDOWN", result[2].replace("-", ""))  # must preserve core

    def test_alignment_with_gaps(self):
        """Test sequences that require gap insertions"""
        seqs = ["ACGT", "AGT", "ACT"]
        aligned = self.feng.align(seqs).splitlines()

        self.assertEqual(len(aligned), 3)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))  # ensure aligned length

        # Confirm original sequences are still in alignment (minus gaps)
        for orig, aligned_seq in zip(seqs, aligned):
            self.assertEqual(orig, aligned_seq.replace("-", ""))

    def test_single_sequence(self):
        """Aligning one sequence returns it unchanged"""
        seq = ["ACTG"]
        result = self.feng.align(seq).splitlines()

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "ACTG")

    def test_empty_sequence_list(self):
        """Handle empty input list"""
        with self.assertRaises(IndexError):
            self.feng.align([])

    def test_supported_algorithms(self):
        """Test that algorithm lists include known methods"""
        pairwise = FengDoolittle.supported_pairwise_algs()
        clustering = FengDoolittle.supported_clustering_algs()

        self.assertIn("needleman_wunsch", pairwise)
        self.assertIn("neighbor_joining", clustering)

    def test_alternate_pairwise(self):
        """Test using a different pairwise alignment algorithm (e.g., gotoh)"""
        gotoh_feng = FengDoolittle(pairwise="gotoh")
        seqs = ["ACGT", "AGT", "ACT"]
        aligned = gotoh_feng.align(seqs).splitlines()

        self.assertEqual(len(aligned), 3)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))

    def test_internal_distance_matrix(self):
        """Check that distances are symmetric and non-negative"""
        seqs = ["ACGT", "AGT", "ACT"]
        _, dist_matrix = self.feng(seqs)

        self.assertEqual(dist_matrix.shape, (3, 3))
        for i in range(3):
            self.assertEqual(dist_matrix[i, i], 0.0)
            for j in range(3):
                self.assertGreaterEqual(dist_matrix[i, j], 0.0)
                self.assertEqual(dist_matrix[i, j], dist_matrix[j, i])

    def test_profile_merge_length(self):
        """Merged profiles should produce equal-length aligned sequences"""
        aligned = self.feng.align(["AAAG", "AG", "AAGG"]).splitlines()

        self.assertTrue(all(len(seq) == len(aligned[0]) for seq in aligned))

    def test_case_insensitivity(self):
        """Ensure alignment behaves case-insensitively (via pairwise alignment)"""
        upper = self.feng.align(["ACTG", "ACTG", "ACTG"])
        lower = self.feng.align(["actg", "actg", "actg"])
        mixed = self.feng.align(["AcTg", "aCtG", "ACTG"])

        self.assertEqual(upper, lower)
        self.assertEqual(lower, mixed)


if __name__ == "__main__":
    unittest.main()
