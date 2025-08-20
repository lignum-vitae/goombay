try:
    # external dependencies
    import numpy
    from numpy import float32
    from numpy._typing import NDArray

    # global packages serving as a placeholder for parsing newick strings - adahik
    from Bio.Phylo.BaseTree import Tree
except ImportError:
    raise ImportError(
        "Please ensure that both numpy and biopython packages are installed.\n"
        "Run the 'pip install <package name>' command from the terminal "
        "if they are not installed or run 'pip install --upgrade <package name>' "
        "to install the latest version."
    )


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
)


from goombay.phylo import NeighborJoining, NewickFormatter

__all__ = ["FengDoolittle", "feng_doolittle"]


def main():
    seq1 = "HOUSEOFCARDSFALLDOWN"
    seq2 = "HOUSECARDFALLDOWN"
    seq3 = "FALLDOWN"
    seq_list = [seq1, seq2, seq3]

    print(feng_doolittle.align([seq1, seq2, seq3]))
    print(FengDoolittle.supported_pairwise_algs())
    print(FengDoolittle.supported_clustering_algs())
    fd_gotoh = FengDoolittle(pairwise="gotoh")
    print(fd_gotoh.align(seq_list))


class FengDoolittle:
    supported_pairwise = {
        "needleman_wunsch": NeedlemanWunsch,
        "jaro": Jaro,
        "jaro_winkler": JaroWinkler,
        "gotoh": Gotoh,
        "wagner_fischer": WagnerFischer,
        "waterman_smith_beyer": WatermanSmithBeyer,
        "hirschberg": Hirschberg,
        "lowrance_wagner": LowranceWagner,
    }

    pw_abbreviations = {
        "nw": "needleman_wunsch",
        "j": "jaro",
        "jw": "jaro_winkler",
        "g": "gotoh",
        "wf": "wagner_fischer",
        "wsb": "waterman_smith_beyer",
        "h": "hirschberg",
        "lw": "lowrance_wagner",
    }

    supported_clustering = {
        "neighbor_joining": NeighborJoining,
    }

    cl_abbreviations = {"nj": "neighbor_joining"}

    def __init__(self, cluster: str = "nj", pairwise: str = "nw", scoring_matrix=None):
        """Initialize Feng-Doolittle algorithm with chosen pairwise method"""
        # Get pairwise alignment algorithm
        if pairwise.lower() in self.supported_pairwise:
            pairwise_class = self.supported_pairwise[pairwise]
        elif pairwise.lower() in self.pw_abbreviations:
            pairwise_class = self.supported_pairwise[self.pw_abbreviations[pairwise]]
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {pairwise}")

        if cluster.lower() in self.supported_clustering:
            self.cluster = self.supported_clustering[cluster]
        elif cluster.lower() in self.cl_abbreviations:
            self.cluster = self.supported_clustering[self.cl_abbreviations[cluster]]
        else:
            raise ValueError(f"Unsupported clustering algorithm: {cluster}")

        if scoring_matrix is not None:
            # getattr(object, attribute, default)
            if getattr(pairwise_class, "supports_scoring_matrix", False):
                self.pairwise = pairwise_class(scoring_matrix=scoring_matrix)
            else:
                raise ValueError(
                    f"The selected pairwise method '{pairwise}' does not support a scoring matrix."
                )
        else:
            self.pairwise = pairwise_class()

    @classmethod
    def supported_pairwise_algs(cls):
        return list(cls.supported_pairwise)

    @classmethod
    def supported_clustering_algs(cls):
        return list(cls.supported_clustering)

    def __call__(self, seqs: list[str]) -> tuple[dict[str, list[str]], NDArray]:
        # tuple[dict[], list[]]:
        """"""
        # This sets the unnormalized sequence distance
        dist_mat_len = len(seqs)
        seq_dist_matrix = numpy.zeros((dist_mat_len, dist_mat_len), dtype=float32)
        profile_dict = {}
        for i, i_seq in enumerate(seqs):
            profile_dict[str(i)] = [i_seq]  # storing lists instead of strings
            for j, j_seq in enumerate(seqs):
                if i < j and i != j:
                    alignment_score = self.pairwise.distance(i_seq, j_seq)
                    seq_dist_matrix[i][j] = alignment_score
                    seq_dist_matrix[j][i] = alignment_score
        return profile_dict, seq_dist_matrix

    def _align(
        self, newick_tree: Tree, profile_dict: dict[str, list[str]], verbose: bool
    ) -> None:
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

    def merge_profiles(self, profile1: list[str], profile2: list[str]) -> list[str]:
        # Pick first seq from each profile as representative (simplified for now)
        rep1 = profile1[0]
        rep2 = profile2[0]

        # Align the two representative sequences
        # Andrew Dahik: Score matrix is here if one wanted to see it.
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


feng_doolittle = FengDoolittle()


if __name__ == "__main__":
    main()
