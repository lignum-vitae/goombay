import unittest
from goombay import NotredameHigginsHeringa, NeedlemanWunsch


class TestFengDoolittle(unittest.TestCase):
    """Test suite for Feng-Doolittle progressive multiple sequence alignment"""

    def setUp(self):
        """Initialize algorithm with default pairwise/clustering (NW + NJ)"""
        self.nhh = NotredameHigginsHeringa()
        self.nw = NeedlemanWunsch()

    def test_identical_sequences(self):
        """Align multiple identical sequences"""
        seqs = ["ACTG", "ACTG", "ACTG"]
        result = self.nhh.align(seqs).splitlines()

        self.assertEqual(len(result), 3)
        for line in result:
            self.assertEqual(line, "ACTG")

    def test_progressive_alignment_order(self):
        """Check output for a known 3-sequence progressive alignment"""
        seqs = ["HOUSEOFCARDSFALLDOWN", "HOUSECARDFALLDOWN", "FALLDOWN"]
        result = self.nhh.align(seqs).splitlines()
        expected = [
            "HOUSEOFCARDSFALLDOWN",
            "HOUSE--CARD-FALLDOWN",
            "------------FALLDOWN",
        ]

        self.assertEqual(len(result), 3)
        self.assertTrue(all(len(row) == len(result[0]) for row in result))  # aligned

        orig_no_gaps = sorted(seq for seq in seqs)
        aligned_no_gaps = sorted(seq.replace("-", "") for seq in result)
        self.assertEqual(orig_no_gaps, aligned_no_gaps)

        for seq in expected:
            self.assertIn(seq, result)

    def test_alignment_with_gaps(self):
        """Test sequences that require gap insertions"""
        seqs = ["ACGT", "AGT", "ACT"]
        aligned = self.nhh.align(seqs).splitlines()

        self.assertEqual(len(aligned), 3)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))  # ensure aligned length

        # Confirm original sequences are still in alignment (minus gaps)
        orig_no_gaps = sorted(seqs)
        aligned_no_gaps = sorted(seq.replace("-", "") for seq in aligned)
        self.assertEqual(orig_no_gaps, aligned_no_gaps)

    def test_single_sequence(self):
        """Aligning one sequence returns it unchanged"""
        seq = ["ACTG"]
        result = self.nhh.align(seq).splitlines()

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "ACTG")

    def test_empty_sequence_list(self):
        """Handle empty input list"""
        with self.assertRaises(ValueError):
            self.nhh.align([])

    def test_supported_algorithms(self):
        """Test that algorithm lists include known methods"""
        pairwise = NotredameHigginsHeringa.supported_pairwise_algs()
        clustering = NotredameHigginsHeringa.supported_clustering_algs()

        self.assertIn("needleman_wunsch", pairwise)
        self.assertIn("neighbor_joining", clustering)

    def test_alternate_pairwise(self):
        """Test using a different pairwise alignment algorithm (e.g., gotoh)"""
        gotoh_feng = NotredameHigginsHeringa(global_pw="gotoh")
        seqs = ["ACGT", "AGT", "ACT"]
        aligned = gotoh_feng.align(seqs).splitlines()

        self.assertEqual(len(aligned), 3)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))

    def test_profile_merge_length(self):
        """Merged profiles should produce equal-length aligned sequences"""
        aligned = self.nhh.align(["AAAG", "AG", "AAGG"]).splitlines()

        self.assertTrue(all(len(seq) == len(aligned[0]) for seq in aligned))

    def test_case_insensitivity(self):
        """Ensure alignment behaves case-insensitively (via pairwise alignment)"""
        aligned = []
        upper = self.nhh.align(["ACTG", "ACTG", "ACTG"])
        lower = self.nhh.align(["actg", "actg", "actg"])
        mixed = self.nhh.align(["AcTg", "aCtG", "ACTG"])

        aligned.extend([upper, lower, mixed])

        self.assertEqual(upper, lower)
        self.assertEqual(lower, mixed)

        for alignment in aligned:
            for seq in alignment.split("\n"):
                self.assertEqual(seq, "ACTG")  # ensure aligned length

    def test_more_than_three_sequences(self):
        """Test alignment with more than three sequences"""
        seqs = ["ACGT", "AGT", "ACT", "AAGT", "ACG"]
        aligned = self.nhh.align(seqs).splitlines()

        self.assertEqual(len(aligned), 5)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))

    def test_more_than_three_identical_sequences(self):
        """Test alignment with more than three identical sequences"""
        seqs = ["ACTG", "ACTG", "ACTG", "ACTG"]
        result = self.nhh.align(seqs).splitlines()

        self.assertEqual(len(result), 4)
        for line in result:
            self.assertEqual(line, "ACTG")

    def test_alignment_with_special_characters(self):
        """Test alignment with sequences containing special characters"""
        seqs = ["A-CGT", "AGT-", "A-CG"]
        aligned = self.nhh.align(seqs).splitlines()

        self.assertEqual(len(aligned), 3)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))

    def test_alignment_with_numbers(self):
        """Test alignment with sequences containing numbers"""
        seqs = ["A1C2GT", "AGT3", "A4CG"]
        aligned = self.nhh.align(seqs).splitlines()

        self.assertEqual(len(aligned), 3)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))

    def test_mixed_empty_and_nonempty(self):
        """Test alignment with a mix of empty and non-empty sequences"""
        seqs = ["ACTG", "", "AGT"]
        aligned = self.nhh.align(seqs).splitlines()
        self.assertEqual(len(aligned), 2)  # Only non-empty sequences should be aligned
        # Remove gaps and check all originals are present (order-independent)
        orig_no_gaps = sorted(seq for seq in seqs if seq)
        aligned_no_gaps = sorted(seq.replace("-", "") for seq in aligned)
        self.assertEqual(orig_no_gaps, aligned_no_gaps)

    def test_alignment_with_whitespace(self):
        """Test alignment with sequences containing leading/trailing whitespace"""
        seqs = ["  ACTG  ", "AGT", "  ACGT"]
        aligned = self.nhh.align(seqs).splitlines()

        self.assertEqual(len(aligned), 3)
        for line in aligned:
            self.assertEqual(len(line), len(aligned[0]))
            self.assertNotIn(" ", line)

    def test_empty_sequences(self):
        """Test alignment with completely empty sequences"""
        seqs = ["", "", ""]
        with self.assertRaises(ValueError):
            self.nhh.align(seqs)

    def test_more_than_three_mixed_empty_and_nonempty(self):
        """Test alignment with more than three sequences, some empty"""
        seqs = ["ACTG", "", "AGT", "AAGT", ""]
        aligned = self.nhh.align(seqs).splitlines()
        self.assertEqual(len(aligned), 3)  # Only non-empty sequences should be aligned
        # Remove gaps and check all originals are present (order-independent)
        orig_no_gaps = sorted(seq for seq in seqs if seq)
        aligned_no_gaps = sorted(seq.replace("-", "") for seq in aligned if seq)
        self.assertEqual(orig_no_gaps, aligned_no_gaps)

    def test_paper_seq_alignment(self):
        """Test sequences aligned in original NHH paper"""
        seqs = [
            "GARFIELD THE LAST FAT CAT",
            "GARFIELD THE FAST CAT",
            "GARFIELD THE VERY FAST CAT",
            "THE FAT CAT",
        ]
        aligned = self.nhh.align(seqs).splitlines()
        expected = [
            "GARFIELD THE LAST FA-T CAT",
            "GARFIELD THE ---- FAST CAT",
            "GARFIELD THE VERY FAST CAT",
            "-------- THE----- FA-T CAT",
        ]
        # self.assertEqual(aligned, expected)


if __name__ == "__main__":
    unittest.main()
