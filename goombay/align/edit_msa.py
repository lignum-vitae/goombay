try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray

    # global packages serving as a placeholder for parsing newick strings - adahik
    from Bio import Phylo, Align
    from io import StringIO
except ImportError:
    raise ImportError("Please pip install all dependencies from requirements.txt!")


# internal dependencies
from goombay.align.edit import (
    NeedlemanWunsch,
    LowranceWagner,
    WagnerFischer,
    WatermanSmithBeyer,
    Gotoh,
    Hirschberg,
    Jaro,
    JaroWinkler,
    SmithWaterman,
    GotohLocal,
)


from goombay.phylo.cluster import (
    NeighborJoining,
    NewickFormatter,
)  # cant run this file directly, use python -m goombay.align.edit_msa in $PWD$/goombay/goombay directory

__all__ = ["FengDoolittle", "feng_doolittle", "NotredameHigginsHeringa", "nhh"]


def main():
    seqs = ["ATCG", "TCG", "ACG"]
    feng_doolittle(seqs)


class MSABase:
    global_supported_pairwise = {
        "needleman_wunsch": NeedlemanWunsch,
        "jaro": Jaro,
        "jaro_winkler": JaroWinkler,
        "gotoh": Gotoh,
        "wagner_fischer": WagnerFischer,
        "waterman_smith_beyer": WatermanSmithBeyer,
        "hirschberg": Hirschberg,
        "lowrance_wagner": LowranceWagner,
    }

    global_pw_abbreviations = {
        "nw": "needleman_wunsch",
        "j": "jaro",
        "jw": "jaro_winkler",
        "gg": "gotoh",
        "wf": "wagner_fischer",
        "wsb": "waterman_smith_beyer",
        "h": "hirschberg",
        "lw": "lowrance_wagner",
    }

    local_supported_pairwise = {
        "smith_waterman": SmithWaterman,
        "gotoh_local": GotohLocal,
    }

    local_pw_abbreviations = {
        "sw": "smith_waterman",
        "gl": "gotoh_local",
    }

    supported_clustering = {
        "neighbor_joining": NeighborJoining,
    }

    cl_abbreviations = {"nj": "neighbor_joining"}

    @classmethod
    def supported_pairwise_algs(cls):
        return list(cls.global_supported_pairwise)

    @classmethod
    def supported_clustering_algs(cls):
        return list(cls.supported_clustering)

    # helper functions for interpreting the newick formatted tree and merging them
    def parse_newick(self, newick: str):
        """
        takes a newick string and converts it into a simple binary tree with Biopythons phylo module

        Args:
        newick: A string in Newick Format containing a MSA that can be inputted into Biopython's Phylo
        """
        tree = Phylo.read(StringIO(newick), "newick")
        return tree

    def merge_profiles(self, profile1: list[str], profile2: list[str]) -> list[str]:
        """
        Docstring for merge_profiles

        :param self: Description
        :param profile1: Description
        :type profile1: list[str]
        :param profile2: Description
        :type profile2: list[str]
        :return: Description
        :rtype: list[str]
        """
        # Pick first seq from each profile as representative (simplified for now)\
        rep1 = profile1[0]
        rep2 = profile2[0]
        # Align the two representative sequences
        aligned_rep1, aligned_rep2 = self.pairwise.align(rep1, rep2).split("\n")

        # Helper: apply gaps to all sequences in a profile
        def apply_gaps(profile, aligned_rep):
            gapped_profile = []
            for seq in profile:
                gapped_seq = []
                seq_i = 0
                for char in aligned_rep:
                    if char == "-":
                        gapped_seq.append("-")
                    else:
                        # Now take the actual character
                        if seq_i < len(seq):
                            gapped_seq.append(seq[seq_i])
                            seq_i += 1

                gapped_profile.append("".join(gapped_seq))
            return gapped_profile

        # Apply alignment gap pattern to all sequences
        aligned_profile1 = apply_gaps(profile1, aligned_rep1)
        aligned_profile2 = apply_gaps(profile2, aligned_rep2)

        return aligned_profile1 + aligned_profile2

    def _align(
        self,
        newick_tree: Phylo.Newick.Tree,
        profile_dict: dict[str, list[str]],
        verbose: bool,
    ):
        """
        Docstring for _align

        :param self: Description
        :param newick_tree: Description
        :param profile_dict: Description
        :param verbose: boolean value for controlling print statements
        :type verbose: bool
        """
        for clade in newick_tree.get_nonterminals(order="postorder"):
            left, right = clade.clades
            if verbose:
                print(left, right)
            if left.name in profile_dict and right.name in profile_dict:
                if verbose:
                    print(f"Merging {left.name} and {right.name}")
                # Merge the aligned profiles
                left_profile = profile_dict[left.name]
                right_profile = profile_dict[right.name]
                merged_profile = self.merge_profiles(
                    left_profile, right_profile
                )  # these are consensus merges

                # Assign unique name to internal node if needed
                if not clade.name:
                    clade.name = f"merged_{id(clade)}"
                # store the merged profile
                profile_dict[clade.name] = merged_profile

    def align(self, seqs: list[str], verbose: bool = False) -> str:
        if not isinstance(seqs, list):
            raise TypeError("Input must be a list of sequences.")

        if not all(isinstance(seq, str) for seq in seqs):
            raise TypeError("All elements in the input list must be strings.")
        seqs = [seq.strip().upper() for seq in seqs if seq.strip()]
        if not seqs:
            raise ValueError("Input list is empty or contains only whitespace strings.")
        if len(seqs) == 1:
            return seqs[0]
        profile_dict, dist_matrix = self(seqs)
        # adding functionality for different clustering algorithms
        nw = self.cluster(dist_matrix).generate_newick()
        newick_tree = NewickFormatter(dist_matrix).parse_newick(nw)
        self._align(newick_tree, profile_dict, verbose)
        if verbose:
            print(dist_matrix)
            print(nw)
            print(newick_tree)

        aligned_seqs = max(profile_dict.values(), key=len)
        rtn_str = []
        for i in range(len(aligned_seqs)):
            rtn_str.append(aligned_seqs[i])
        return "\n".join(rtn_str)

    def get_matrix(self, seqs: list[str], verbose: bool = False):  # getter class
        return self(seqs)[1]  # should return distance matrix

    def gen_profile_dict(self, seqs: list[str]):
        profile_dict = {}
        for i in range(len(seqs)):
            profile_dict[str(i)] = [seqs[i]]  # storing lists instead of strings
        return profile_dict

    def _track_sequences(self, seqs: list[str]):
        """
        This is a utility function for tracking the order of the sequences as they were inputed from their dataset (list)

        :param self: Description
        :param seqs: list(str) -> a list of strings that represent sequences that are retrieved from their respective index

        returns: a library with key-value pairs -> [key: chr(seqs.index(i_seq) + 65).lower())] value: i_seq
        """
        lib_seq = {}
        for i in range(len(seqs)):
            lib_seq[str(i)] = seqs[i]
        return lib_seq

    def _create_positions(self, seq_x):
        seq_x_positions = []
        count = 0
        for char in seq_x:
            if char != "-":
                count += 1
                seq_x_positions.append(count)
            else:
                seq_x_positions.append("-")
        return seq_x_positions

    def _rev_pair(self, pair):
        return (pair[-1], pair[0])


class FengDoolittle(MSABase):
    """functions below are unique to FengDoolittle"""

    def __init__(self, cluster: str = "nj", pairwise: str = "nw"):
        """Initialize Feng-Doolittle algorithm with chosen pairwise method"""
        # Get pairwise alignment algorithm
        if pairwise.lower() in self.global_supported_pairwise:
            self.pairwise = self.global_supported_pairwise[pairwise]()
        elif pairwise.lower() in self.global_pw_abbreviations:
            self.pairwise = self.global_supported_pairwise[
                self.global_pw_abbreviations[pairwise]
            ]()
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {pairwise}")

        if cluster.lower() in self.supported_clustering:
            self.cluster = self.supported_clustering[cluster]
        elif cluster.lower() in self.cl_abbreviations:
            self.cluster = self.supported_clustering[self.cl_abbreviations[cluster]]
        else:
            raise ValueError(f"Unsupported clustering algorithm: {cluster}")

    def __call__(self, seqs: list[str]):
        """"""
        # This sets the unnormalized sequence distance
        dist_mat_len = len(seqs)
        seq_dist_matrix = numpy.zeros((dist_mat_len, dist_mat_len), dtype=float64)

        # FengDoolittle can make use of a different algorithm to generate both the profile_dict and sequence distance matrix
        # Might be worth moving it to a function

        # storing lists instead of strings
        profile_dict = {str(i): [seq] for i, seq in enumerate(seqs)}
        for i, i_seq in enumerate(seqs):
            for j, j_seq in enumerate(seqs):
                if i < j and i != j:
                    alignment_score = self.pairwise.distance(i_seq, j_seq)
                    seq_dist_matrix[i][j] = alignment_score
                    seq_dist_matrix[j][i] = alignment_score
        return profile_dict, seq_dist_matrix

    def align(self, seqs: list[str], verbose: bool = False) -> str:
        return super().align(seqs, verbose)


class NotredameHigginsHeringa(MSABase):  # T-Coffee implementation
    """functions below are unique to NotredameHigginsHeringa"""

    # instead of waterman smith bayer and needleman wunsch, follow RNAinformatik's use of Gotoh
    def __init__(
        self, local_pw: str = "sw", global_pw: str = "nw", cluster: str = "nj"
    ):
        """Initialize NotredameHigginsHeringa algorithm with chosen methods"""
        # Get Global pairwise alignment algorithm
        if global_pw.lower() in self.global_supported_pairwise:
            self.global_pw = self.global_supported_pairwise[global_pw]()
        elif global_pw.lower() in self.global_pw_abbreviations:
            self.global_pw = self.global_supported_pairwise[
                self.global_pw_abbreviations[global_pw]
            ]()
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {global_pw}")
        # get local pairwise alignment algorithm
        if local_pw.lower() in self.local_supported_pairwise:
            self.local_pw = self.local_supported_pairwise[local_pw]()
        elif local_pw.lower() in self.local_pw_abbreviations:
            self.local_pw = self.local_supported_pairwise[
                self.local_pw_abbreviations[local_pw]
            ]()
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {local_pw}")

        # self.pairwise is set equal to the global_pw string
        self.pairwise = self.global_pw

        # get clustering algorithm - one the extended library use a clustering method to make a di
        if cluster.lower() in self.supported_clustering:
            self.cluster = self.global_supported_pairwise[cluster]()
        elif cluster.lower() in self.cl_abbreviations:
            self.cluster = self.supported_clustering[self.cl_abbreviations[cluster]]
        else:
            raise ValueError(f"Unsupported clustering algorithm: {cluster}")

    # Andrew Dahik: we may want to move the primary library into its own file like edit distance!
    # heuristic consistency score
    def __call__(self, seqs: list[str]):
        """
        Instead of distance, using heuristic consistency score
        Args:
        seqs (list[str]) - an input of given sequences for alignment

        Returns:
        Profile Dictionary or history of sequences and Distance Matrix
        """
        # Similar to FengDoolittle, Store a profile dict to leverage align
        profile_dict = self.gen_profile_dict(seqs)
        # first perform optimum pairwise alignments for each sequence against each other
        seq_tracker = self._track_sequences(seqs)
        alignments = self.compute_alignments(seqs)

        primary_library = self.compute_primary_library(alignments, seqs, seq_tracker)

        """Extension library logic"""
        # taking the primary library
        # for each stored sequence need to calculate extended library score against other alignments

        # alrighty, lets do the extension
        extended_library = self.form_extension_library(primary_library, len(seqs))

        # distance matrix from the extended library, almost ready to be clustered
        seq_dist_matrix = self.create_distance_matrix_NNH(extended_library, seq_tracker)
        return profile_dict, seq_dist_matrix

    def compute_alignments(self, seqs: list[str]):
        r"""
        Docstring for merge_alignments: this is a pairwise optimal alignment step following the instruction from https://backofenlab.github.io/BioinformaticsII-pages/exercise-sheet-3.html

        Args:
        :param self: for object instantiation calling
        :param seqs: list(str) -> this is the input of strings to be aligned such that A(x,y) for all sequences a,...,n

        Returns:
        a library of merged alignments with the following key-value format:
         {'a-b': (sequence_a_aligned\sequence_b_aligned, distance score}
        the key is a string that shows the ith sequence given with its corresponding unicode character (the first sequence is a, second b, ...) aligned (denoted with -) with the comparing sequence
        """
        # maybe add a counter for easier to track keys?
        # This sets the unnormalized sequence distance, odd numbered keys are aligned to global alignments, even are assigned to local alignments
        merged_alignments = {}
        aligner = Align.PairwiseAligner(match_score=1.0)
        for i, i_seq in enumerate(seqs):
            for j, j_seq in enumerate(seqs):
                # this condition helps avoid repeats
                if i < j and i != j:
                    """"""
                    # use Biopythons pairwise align for now to compare alignments
                    alignment = aligner.align(i_seq, j_seq)
                    aligned_seq_ij = alignment[0][0] + "\n" + alignment[0][1]

                    # DEBUG: Check alignment lengths
                    lines = aligned_seq_ij.split("\n")
                    seq1_aligned = lines[0]
                    seq2_aligned = lines[1]

                    key = f"{i}-{j}"
                    merged_alignments[key] = (
                        aligned_seq_ij,
                        self._merge_weight_calc(aligned_seq_ij),
                    )
                    # local_align = self.local_pw.align(i_seq, j_seq)
                    # if global_align != local_align:
                    # merged_alignments[
                    # chr(seqs.index(i_seq) + 65).lower()
                    # + "-"
                    # + chr(seqs.index(j_seq) + 65).lower()
                    # ] = (local_align, self._merge_weight_calc(local_align))
        return merged_alignments

    def _merge_weight_calc(self, seq_ab):
        """
        Docstring for _merge_weight_calc

        Args:
        :param self: for object instantiation calling
        :param seq_ab: seq_ab is the A(a,b) or alignment of sequence a and sequence b
        """
        total_match = 0
        sequences = seq_ab.split("\n")
        seq_a, seq_b = sequences[0], sequences[1]
        # get the total character length not including gaps of the min seq
        seqa_trimmed, seqb_trimmed = seq_a.replace("-", ""), seq_b.replace("-", "")
        total_count = min(len(seqa_trimmed), len(seqb_trimmed))
        for i in range(len(seq_a)):
            if seq_a[i] != "-" and seq_b[i] != "-":  # skip over indels if both are gaps
                if seq_a[i] == seq_b[i]:
                    total_match += 1
        return round((total_match * 100) / (total_count), 1)

    def compute_primary_library(
        self,
        alignments: dict[str, tuple[str, float]],
        seqs: list[str],
        seq_tracker: dict[str],
    ):
        # primary library forming logic - primary library is interesting in that its actually a dictionary containing a list of tuples as the values assigned to keys determined by their indexes unicode value
        primary_library = {}
        for m in range(len(seqs)):
            for n in range(m + 1, len(seqs)):
                if n < len(seqs):
                    primary_library[chr(m + 65).lower() + "-" + chr(n + 65).lower()] = (
                        self._form_primary_constraints(alignments, seq_tracker, m, n)
                    )
        return primary_library

    def _form_primary_constraints(
        self,
        alignments: dict[str, tuple[str, float]],
        seq_tracker: dict[str, str],
        m: int,
        n: int,
    ):  # forming the primary library
        """
        Docstring for form_primary_constraints: follows the examples from https://backofenlab.github.io/BioinformaticsII-pages/exercise-sheet-3.html

        :param self: Description
        :param seqs: Original input, to hold the original positions when forming the library
        :param A: Alignments, calculated from compute_alignments

        """
        seq_keys, alignment_keys = list(seq_tracker.keys()), list(alignments.keys())
        primary_lib_tuples = []

        # The primary sequence is the first one
        seqi_key = seq_keys[m]
        orig_seqi = seq_tracker[seqi_key][:]

        # Compare against sequence 'b' (second sequence)
        seqj_key = seq_keys[n]
        orig_seqj = seq_tracker[seqj_key][:]

        # Get the alignment between primary and secondary
        alignment_key = f"{seqi_key}-{seqj_key}"
        current_aligned_seqs = alignments[alignment_key][0].split("\n")
        aligned_primary = current_aligned_seqs[0]
        aligned_secondary = current_aligned_seqs[1]
        # get the weight for each compared position
        weight = alignments[alignment_key][1]

        # Track positions in original sequences
        i_pos = 0  # position in original primary sequence (0-indexed)
        j_pos = 0  # position in original secondary sequence (0-indexed)

        # Create mapping: which primary positions align to which secondary positions
        alignment_map = {}  # key: (i_pos+1, j_pos+1), value: weight

        for k in range(len(aligned_primary)):
            # booleans for tracking true characters from the original unaligned sequence imposed with its aligned counterpart
            has_primary = aligned_primary[k] != "-"
            has_secondary = aligned_secondary[k] != "-"

            # If both positions have characters and they match, record the alignment
            if has_primary and has_secondary:
                if aligned_primary[k] == aligned_secondary[k]:
                    alignment_map[(i_pos + 1, j_pos + 1)] = weight

            # Increment position counters (only when not a gap)
            if has_primary:
                i_pos += 1
            if has_secondary:
                j_pos += 1

        # Now generate tuples for ALL position combinations in order
        for i in range(1, len(orig_seqi) + 1):
            for j in range(1, len(orig_seqj) + 1):
                weight = alignment_map.get((i, j), 0)
                primary_lib_tuples.append((i, j, weight))

        return primary_lib_tuples

    def form_extension_library(
        self, prim_lib: dict[list[tuple]], n: int
    ):  # n is the number of sequences
        """
        Docstring for form_extension_library
        This method performs the step to increment
        :param self: Description
        :param prim_lib: Description
        :param n: Description
        """
        extension_library = {pair: {} for pair in prim_lib}
        # for tracking iterations, the letters are indexed back to their chr value

        # ab_key represents sequence a and sequence b
        for ab_key in prim_lib.keys():
            a, b = ab_key.split("-")
            possible_pairs = {}
            # initialize all possible pairs identified in primary libary
            for i, j, w_ab in prim_lib[ab_key]:
                possible_pairs[i, j] = w_ab
            # iterating over the tuples stored, which is residue i and residue j and weight from the key
            # find triangulated matchs
            for c in range(n):
                if chr(c + ord("a")) in (a, b):
                    continue  # exit this iteration

                # triangulation step: check all positions k of c, i2 is the same residue as i
                c_key = chr(c + ord("a"))
                ac_key = f"{a}-{c_key}"
                cb_key = f"{c_key}-{b}"
                ca_key = f"{c_key}-{a}"
                bc_key = f"{b}-{c_key}"

                # Get matches from primary library (check both orderings)
                ac_matches = prim_lib.get(ac_key, [])
                ca_matches = [
                    (k, i, w) for i, k, w in prim_lib.get(ca_key, [])
                ]  # Swap i,k for reversed pair

                cb_matches = prim_lib.get(cb_key, [])
                bc_matches = [
                    (j, k, w) for k, j, w in prim_lib.get(bc_key, [])
                ]  # Swap k,j for reversed pair

                # Combine both directions
                all_ac = ac_matches + ca_matches
                all_cb = cb_matches + bc_matches

                # For all matches a[i] -> c[k]
                for i, k, w_ac in all_ac:
                    # For all matches c[k] -> b[j]
                    for k2, j, w_cb in all_cb:
                        if k2 == k:  # Same position k in sequence c
                            # Triangulation found: a[i] -> c[k] -> b[j]
                            triangulation_weight = min(w_ac, w_cb)

                            # Add to the (i,j) pair's total weight
                            if (i, j) not in possible_pairs:
                                possible_pairs[(i, j)] = 0  # No primary library match
                            possible_pairs[(i, j)] += triangulation_weight

            # Store all pairs with their accumulated weights
            extension_library[ab_key] = possible_pairs
        return extension_library

    def create_distance_matrix_NNH(
        self, extension_library: dict[dict], seq_tracker: dict[str]
    ):
        """
        Docstring for create_distance_matrix_NNH

        :param self: Description
        :param extension_library: Description
        :type extension_library: dict[dict]
        :param seqs: Description
        :type seqs: list[str]
        """
        n = len(seq_tracker)
        distance_matrix = numpy.zeros((n, n))

        for i, a in enumerate(seq_tracker.keys()):
            for j, b in enumerate(seq_tracker.keys()):
                if i >= j:  # Skip diagonal and lower triangle
                    continue

                # Try both orderings of the pair
                pair_key = f"{a}-{b}"
                reverse_key = f"{b}-{a}"

                if pair_key in extension_library:
                    pair_dict = extension_library[pair_key]
                elif reverse_key in extension_library:
                    pair_dict = extension_library[reverse_key]
                else:
                    continue

                # Collect all non-zero scores
                scores = [
                    weight for (ai, bj), weight in pair_dict.items() if weight > 0
                ]

                if not scores:
                    distance_matrix[i, j] = 1.0
                    distance_matrix[j, i] = 1.0
                else:
                    avg_score = numpy.mean(scores)
                    similarity = avg_score / 100.0
                    distance = 1.0 - similarity

                    distance_matrix[i, j] = distance
                    distance_matrix[j, i] = distance
        return distance_matrix


feng_doolittle = FengDoolittle()
nhh = NotredameHigginsHeringa()

if __name__ == "__main__":
    main()
