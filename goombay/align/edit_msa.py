try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray

    # global packages serving as a placeholder for parsing newick strings - adahik
    from Bio import Phylo
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
    GotohLocal,
)


from goombay.phylo.cluster import NeighborJoining, NewickFormatter

__all__ = ["FengDoolittle", "feng_doolittle", "TCoffee", "t_coffee"]


def main():
    # FengDoolittle
    # seq1 = "HOUSEOFCARDSFALLDOWN"
    # seq2 = "HOUSECARDFALLDOWN"
    # seq3 = "FALLDOWN"
    # seq_list = [seq1, seq2, seq3]

    # print(FengDoolittle.align(seq_list))
    # print(FengDoolittle.supported_pairwise_algs())
    # print(FengDoolittle.supported_clustering_algs())
    # fd_gotoh = FengDoolittle(pairwise="gotoh")

    # TCoffee seqs from rna.informatik
    seq1 = "GARFIELD_THE_LAST_FAT_CAT"
    seq2 = "GARFIELD_THE_FAST_CAT"
    seq3 = "GARFIELD_THE_VERY_FAST_CAT"
    seq4 = "THE_FAT_CAT"
    seq_list = [seq1, seq2, seq3, seq4]
    print(TCoffee.supported_pairwise_algs())
    print(TCoffee.supported_clustering_algs())
    print(t_coffee(seq_list))
    # fd_gotoh = TCoffee(pairwise="gotoh")
    # print(t_coffee.align(seq_list))


class AlignTemp:
    supported_pairwise = {
        "needleman_wunsch": NeedlemanWunsch,
        "jaro": Jaro,
        "jaro_winkler": JaroWinkler,
        "gotoh": Gotoh,
        "wagner_fischer": WagnerFischer,
        "waterman_smith_beyer": WatermanSmithBeyer,
        "hirschberg": Hirschberg,
        "lowrance_wagner": LowranceWagner,
        "gotoh_local": GotohLocal,
    }

    pw_abbreviations = {
        "nw": "needleman_wunsch",  # global alignment
        "j": "jaro",
        "jw": "jaro_winkler",
        "gg": "gotoh",
        "wf": "wagner_fischer",
        "wsb": "waterman_smith_beyer",  # local alignment
        "h": "hirschberg",
        "lw": "lowrance_wagner",
        "gl": "gotoh_local",
    }

    supported_clustering = {
        "neighbor_joining": NeighborJoining,
    }

    cl_abbreviations = {"nj": "neighbor_joining"}

    @classmethod
    def supported_pairwise_algs(cls):
        return list(cls.supported_pairwise)

    @classmethod
    def supported_clustering_algs(cls):
        return list(cls.supported_clustering)

    # helper functions for interpreting the newick formatted tree and merging them
    def parse_newick(self, newick):
        """takes a newick string and converts it into a simple binary tree with Biopythons phylo module"""
        tree = Phylo.read(StringIO(newick), "newick")
        return tree

    def merge_profiles(self, profile1: list[str], profile2: list[str]) -> list[str]:
        # Pick first seq from each profile as representative (simplified for now)
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
                        gapped_seq.append(seq[seq_i])
                        seq_i += 1
                gapped_profile.append("".join(gapped_seq))
            return gapped_profile

        # Apply alignment gap pattern to all sequences
        aligned_profile1 = apply_gaps(profile1, aligned_rep1)
        aligned_profile2 = apply_gaps(profile2, aligned_rep2)

        return aligned_profile1 + aligned_profile2

    def _align(self, newick_tree, profile_dict, verbose: bool):
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


class FengDoolittle(AlignTemp):
    """functions beloww are unique to FengDoolittle"""

    def __init__(self, cluster: str = "nj", pairwise: str = "nw"):
        """Initialize Feng-Doolittle algorithm with chosen pairwise method"""
        # Get pairwise alignment algorithm
        if pairwise.lower() in self.supported_pairwise:
            self.pairwise = self.supported_pairwise[pairwise]()
        elif pairwise.lower() in self.pw_abbreviations:
            self.pairwise = self.supported_pairwise[self.pw_abbreviations[pairwise]]()
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
        profile_dict = {}
        for i, i_seq in enumerate(seqs):
            profile_dict[str(i)] = [i_seq]  # storing lists instead of strings
            for j, j_seq in enumerate(seqs):
                if i < j and i != j:
                    alignment_score = self.pairwise.distance(i_seq, j_seq)
                    seq_dist_matrix[i][j] = alignment_score
                    seq_dist_matrix[j][i] = alignment_score
        return profile_dict, seq_dist_matrix


class TCoffee(AlignTemp):  # Notredame-Higgins-Heringa implementation
    """functions beloww are unique to T-Coffee"""

    # instead of waterman smith bayer and needleman wunsch, follow RNAinformatik's use of Gotoh
    def __init__(
        self, local_pw: str = "wsb", global_pw: str = "nw", cluster: str = "nj"
    ):
        """Initialize T-Coffee algorithm with chosen methods"""
        # Get Global pairwise alignment algorithm
        if global_pw.lower() in self.supported_pairwise:
            self.global_pw = self.supported_pairwise[global_pw]()
        elif global_pw.lower() in self.pw_abbreviations:
            self.global_pw = self.supported_pairwise[self.pw_abbreviations[global_pw]]()
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {global_pw}")
        # get local pairwise alignment algorithm
        if local_pw.lower() in self.supported_pairwise:
            self.local_pw = self.supported_pairwise[local_pw]()
        elif local_pw.lower() in self.pw_abbreviations:
            self.local_pw = self.supported_pairwise[self.pw_abbreviations[local_pw]]()
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {local_pw}")

        # get clustering algorithm - one the extended library use a clustering method to make a di
        if cluster.lower() in self.supported_clustering:
            self.cluster = self.supported_pairwise[cluster]()
        elif cluster.lower() in self.cl_abbreviations:
            self.cluster = self.supported_clustering[self.cl_abbreviations[cluster]]
        else:
            raise ValueError(f"Unsupported clustering algorithm: {cluster}")

    # Andrew Dahik: we may want to move the primary library into its own file like edit distance!
    # heuristic consistency score
    def __call__(self, seqs: list[str]):
        """Instead of distance, using heuristic consistency score"""

        # make merge alignments for tracking the alignments themselves
        # the index positions should keep track
        # REMOVE Print Statements when done
        # print(seqs)
        # print(self._sum_of_pairs(seqs))

        prelim_libr = self.merge_alignments(seqs)

        # print(prelim_libr)
        # placeholder - based on the weighted scores of each alignment individually, hold the best aligned sequence based on the

        best_seqs = self.store_best_sequences(prelim_libr)
        print(best_seqs)
        # Similar to FengDoolittle, Store a profile dict to leverage align
        profile_dict = {}
        for i, i_seq in enumerate(best_seqs):
            profile_dict[str(i)] = [i_seq]  # storing lists instead of strings
        # primary library forming logic - primary library is interesting in that its actually a dictionary containing a list of tuples as the values assigned to keys determined by their indexes unicode value
        primary_library = self.form_primary_constraints(prelim_libr)
        # print(primary_library)
        """Extension library logic"""
        # taking the primary library
        # for each stored sequence need to calculate extended library score against other alignments

        # alrighty, lets do the extension
        extended_library = self.form_extension_library(primary_library, len(seqs))
        # print(extended_library)

        # distance matrix from the extended library, almost ready to be clustered
        seq_dist_matrix = self.create_distance_matrix(
            extended_library, list(best_seqs.keys())
        )
        # print(dist_matrix)
        return profile_dict, seq_dist_matrix

    def merge_alignments(self, seqs):
        # maybe add a counter for easier to track keys?
        # This sets the unnormalized sequence distance, odd numbered keys are aligned to global alignments, even are assigned to local alignments
        merged_alignments = {}
        for i, i_seq in enumerate(seqs):
            for j, j_seq in enumerate(seqs):
                # this condition helps avoid repeats
                if i < j and i != j:
                    global_align = self.global_pw.align(i_seq, j_seq)
                    merged_alignments[
                        chr(seqs.index(i_seq) + 65).lower()
                        + "-"
                        + chr(seqs.index(j_seq) + 65).lower()
                    ] = (global_align, self._merge_weight_calc(global_align))

                    local_align = self.local_pw.align(i_seq, j_seq)
                    if global_align != local_align:
                        merged_alignments[
                            chr(seqs.index(i_seq) + 65).lower()
                            + "-"
                            + chr(seqs.index(j_seq) + 65).lower()
                        ] = (local_align, self._merge_weight_calc(local_align))
        # print(merged_alignments)
        return merged_alignments

    def _merge_weight_calc(self, seq_ab):
        # seq_ab is the A(a,b) - alignment of sequence a and sequence b
        # print("weight calc")
        total_match = 0
        total_mismatch = 0
        sequences = seq_ab.split("\n")
        seq_a, seq_b = sequences[0], sequences[1]
        for i in range(len(seq_a)):
            if seq_a[i] != "-" and seq_b[i] != "-":  # skip over indels
                if seq_a[i] == seq_b[i]:
                    total_match += 1
                else:
                    total_mismatch += 1

        return round((total_match * 100) / (total_match + total_mismatch), 1)

    def store_best_sequences(self, prelim_libr):
        """
        This function works to take the best alignment for each sequence from the original input base on weight
        The input structure of prelim_libr is {a-b: alignment_seq_a\nalignment_seq_b, ... }
        """
        best_sequences = []
        seq_weight_map = {}
        # iterate through the map and store the alignments in a new arrangement
        for ab_key, (align_seq, w_ab) in prelim_libr.items():
            a, b = ab_key.split("-")
            if a not in seq_weight_map:
                seq_weight_map[a] = {}
            if b not in seq_weight_map:
                seq_weight_map[b] = {}
            seq_a, seq_b = align_seq.split("\n")
            if seq_a not in seq_weight_map[a]:
                seq_weight_map[a][seq_a] = w_ab
            else:
                if w_ab > seq_weight_map[a][seq_a]:
                    seq_weight_map[a][seq_a] = w_ab
            if seq_b not in seq_weight_map[a]:
                seq_weight_map[b][seq_b] = w_ab
            else:
                if w_ab > seq_weight_map[b][seq_b]:
                    seq_weight_map[b][seq_b] = w_ab
        # now chose the best aligned sequence for each sequence based on weight value, iterating over sequence character index
        for seq_chr_ind in seq_weight_map:
            # appending to best_sequences a placeholder value that is replaced by the seq value below
            best_sequences.append(0)
            best_score, candidates = 0, seq_weight_map[seq_chr_ind]
            for seq in candidates:
                if candidates[seq] > best_score:
                    best_score = candidates[seq]
                    best_sequences[-1] = seq
        # now take indices again, I'm getting a little lazy here and wiping seq_map since we don't need it any more at this line
        seq_weight_map = {}
        for i in range(len(best_sequences)):
            i_char = chr(i + 65).lower()  # index character i
            seq_weight_map[i_char] = best_sequences[i]
        return seq_weight_map

    def form_primary_constraints(self, prelim_libr):  # forming the primary library
        primary_lib_tuples = []

        for key in prelim_libr:
            seq_ab, weight = prelim_libr[key][0], prelim_libr[key][1]
            seq_a, seq_b = seq_ab.split("\n")

            pos_a_lst, pos_b_lst = self._create_positions(
                seq_a
            ), self._create_positions(seq_b)
            best_matches = {}  # key: a_pos_index, value: (b_pos_index, weight)

            for i in range(len(seq_a)):
                if seq_a[i] == "-":
                    continue

                for j in range(len(seq_b)):
                    if seq_b[j] == "-":
                        continue
                    if seq_a[i] == seq_b[j]:
                        w = 0.0 if seq_a[i] == "-" else weight

                        # pick the best match: higher weight or closer position if weights equal
                        if (
                            i not in best_matches
                            or w > best_matches[i][1]
                            or (
                                w == best_matches[i][1]
                                and abs(i - j) < abs(i - best_matches[i][0])
                            )
                        ):
                            best_matches[i] = (j, w)

            # Append only the best matches
            for i, (j, w) in best_matches.items():
                primary_lib_tuples.append(
                    (
                        key[: key.index("-")],
                        pos_a_lst[i],
                        key[key.index("-") + 1 :],
                        pos_b_lst[j],
                        w,
                    )
                )

        # Compress tuples into final sequences_maps
        sequences_maps = {}
        seq_pair_mem = []
        for tup in primary_lib_tuples:
            seq_pair = (tup[0], tup[2])
            if (
                seq_pair not in seq_pair_mem
                and self._rev_pair(seq_pair) not in seq_pair_mem
            ):
                seq_pair_mem.append(seq_pair)
                sequences_maps[seq_pair[0] + "-" + seq_pair[1]] = []
            sequences_maps[seq_pair[0] + "-" + seq_pair[1]].append(
                (tup[1], tup[3], tup[4])
            )

        # Trim triangulations not starting at pos_b = 1
        for seq_pair in sequences_maps:
            seq_pair_lst = sequences_maps[seq_pair]
            tup_idx = 0
            for tup in seq_pair_lst:
                if tup[1] == 1:
                    break
                tup_idx += 1
            sequences_maps[seq_pair] = seq_pair_lst[tup_idx:]

        return sequences_maps

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

    def form_extension_library(self, prim_lib, n):  # n is the number of sequences
        extension_library = {pair: {} for pair in prim_lib}
        # for tracking iterations, the letters are indexed back to their chr value

        # ab_key represents sequence a and sequence b
        for ab_key, ab_residue_matches in prim_lib.items():
            a, b = ab_key.split("-")
            # iterating over the tuples stored, which is residue i and residue j and weight from the key
            for i, j, w_ab in ab_residue_matches:
                for c in range(n):
                    total_weight = w_ab
                    if chr(c + 65) in (a, b):
                        continue  # exit this iteration
                    # check all positions k of c, i2 is the same residue as i
                    for i2, k, w_ac in prim_lib.get(f"{a}-{c}", []):
                        if i2 == i:
                            for k2, j2, w_cb in prim_lib.get(f"{c}-{b}", []):
                                if (
                                    k2 == k and j2 == j
                                ):  # triangulation: residue at position k is the same for a-c and b-c, since position i and j have the same value, i can associate with k like j can with k
                                    total_weight += min(w_ac, w_cb)
                    # store the final score
                    extension_library[ab_key].setdefault(
                        (i, j), total_weight
                    )  # using pythons dictionary forming method, adding keys
        return extension_library

    def create_distance_matrix(self, extension_library, seq_names):  # TODO
        """
        Taking the output from store_best_sequences as seq_names (a, b, c, ...) and the extension library these seq names are tracked too
        """
        n = len(seq_names)
        # time to make the distance matrix
        distance_matrix = numpy.zeros((n, n))
        for i, a in enumerate(seq_names):
            for j, b in enumerate(seq_names):
                if i == j or a + "-" + b not in extension_library:
                    continue
                scores = []
                # the keys are tuples with the index positions where aligned
                for ai, bj in extension_library[a + "-" + b]:
                    scores.append(extension_library[a + "-" + b][(ai, bj)])
                # average (or sum ) the scores
                avg_scores = numpy.mean(scores) if scores else 0
                # convert similarity to distance
                distance_matrix[i, j] = 1 - avg_scores
                distance_matrix[j, i] = 1 - avg_scores

        return distance_matrix


feng_doolittle = FengDoolittle()
t_coffee = TCoffee()
if __name__ == "__main__":
    main()
