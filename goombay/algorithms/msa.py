# internal dependencies
from goombay.algorithms.editdistance import (
    NeedlemanWunsch,
    LowranceWagner,
    WagnerFischer,
    WatermanSmithBeyer,
    Gotoh,
    Hirschberg,
    Jaro,
    JaroWinkler,
)

try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
except ImportError:
    raise ImportError("Please pip install all dependencies from requirements.txt!")


def main():
    seq1 = "HOUSEOFCARDSFALLDOWN"
    seq2 = "HOUSECARDFALLDOWN"
    seq3 = "FALLDOWN"

    print(longest_common_substring_msa.align([seq2, seq1, seq3]))
    print(longest_common_substring_msa.distance([seq2, seq1, seq3]))
    print(longest_common_substring_msa.normalized_distance([seq2, seq1, seq3]))
    print(longest_common_substring_msa.similarity([seq2, seq1, seq3]))
    print(longest_common_substring_msa.normalized_similarity([seq2, seq1, seq3]))


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

    abbreviations = {
        "nw": "needleman_wunsch",
        "j": "jaro",
        "jw": "jaro_winkler",
        "g": "gotoh",
        "wf": "wagner_fischer",
        "wsb": "waterman_smith_beyer",
        "h": "hirschberg",
        "lw": "lowrance_wagner",
    }

    def __init__(self, pairwise: str = "needleman_wunsch"):
        """Initialize Feng-Doolittle algorithm with chosen pairwise method"""
        # Get pairwise alignment algorithm
        if pairwise in self.supported_pairwise:
            self.pairwise = self.supported_pairwise[pairwise]()
        elif pairwise in self.abbreviations:
            self.pairwise = self.supported_pairwise[self.abbreviations[pairwise]]()
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {pairwise}")

    @classmethod
    def supported_pairwise_algs(cls):
        return list(cls.supported_pairwise)

    def __call__(self):
        raise NotImplementedError("Class method not yet implemented")

    def matrix(self):
        raise NotImplementedError("Class method not yet implemented")

    def align(self):
        raise NotImplementedError("Class method not yet implemented")

    def distance(self):
        raise NotImplementedError("Class method not yet implemented")

    def similarity(self):
        raise NotImplementedError("Class method not yet implemented")

    def normalized_distance(self):
        raise NotImplementedError("Class method not yet implemented")

    def normalized_similarity(self):
        raise NotImplementedError("Class method not yet implemented")


class LongestCommonSubstringMSA:
    def _common_substrings(self, seqs: list[str]) -> list[str]:
        if not isinstance(seqs, list):
            raise TypeError("common_substrings expects a list of strings")
        if len(seqs) != 2:
            raise ValueError("common_substrings requires exactly two strings")

        s1, s2 = seqs
        common = set()
        s2_length = len(s2)

        for start in range(s2_length):
            for end in range(start + 2, s2_length + 1):  # substrings of length >= 2
                substr = s2[start:end]
                if substr not in s1:
                    break
                common.add(substr)

        return list(common)

    def __call__(self, seqs: list[str]) -> list[str]:
        if (not isinstance(seqs, list) and not isinstance(seqs, tuple)) or not all(
            isinstance(s, str) for s in seqs
        ):
            raise TypeError("longest_common_substring_msa expects a list of strings")
        if len(seqs) < 2:
            raise ValueError("Provide at least two seqs")

        seqs = [seq.upper() for seq in seqs]

        # Generate substrings from first and last strings
        motifs = self._common_substrings([seqs[0], seqs[-1]])
        if not motifs:
            return [""]

        longest = []
        longest_len = -1
        motifs.sort(key=len, reverse=True)
        for motif in motifs:
            motif_len = len(motif)
            if all(motif in seq for seq in seqs) and motif_len >= longest_len:
                longest.append(motif)
                longest_len = motif_len
            if motif_len < longest_len:
                break
        return longest

    def align(self, seqs: list[str]) -> list[str]:
        return self(seqs)

    def distance(self, seqs: list[str]) -> int:
        max_match = len(max(seqs, key=len))
        if any(not seq for seq in seqs):
            return max_match
        if max_match <= 1 and all(seqs[0] in seq for seq in seqs):
            return 0
        return max_match - self.similarity(seqs)

    def similarity(self, seqs: list[str]) -> int:
        lcsub = self(seqs)
        if not lcsub:
            return 0
        if all(len(seq) == 0 for seq in seqs) or all(seqs[0] == seq for seq in seqs):
            return len(seqs[0]) if len(seqs[0]) >= 1 else 1
        return len(lcsub[0])

    def normalized_distance(self, seqs: list[str]) -> float:
        return 1 - self.normalized_similarity(seqs)

    def normalized_similarity(self, seqs: list[str]) -> float:
        if len(max(seqs, key=len)) == 1 and all(seqs[0] in seq for seq in seqs):
            return 1.0
        lcsub = self(seqs)
        if not lcsub:
            return 0.0
        return len(lcsub[0]) / len(min(seqs, key=len))


feng_doolittle = FengDoolittle()
longest_common_substring_msa = LongestCommonSubstringMSA()


if __name__ == "__main__":
    main()
