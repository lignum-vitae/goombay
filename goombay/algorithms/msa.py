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

    print(longest_common_substring.align([seq2, seq1, seq3]))
    print(longest_common_substring.distance([seq2, seq1, seq3]))
    print(longest_common_substring.normalized_distance([seq2, seq1, seq3]))
    print(longest_common_substring.similarity([seq2, seq1, seq3]))
    print(longest_common_substring.normalized_similarity([seq2, seq1, seq3]))


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
    def _common_substrings(self, sequences: list[str]) -> list[str]:
        if not isinstance(sequences, list):
            raise TypeError("common_substrings expects a list of strings")
        if len(sequences) != 2:
            raise ValueError("common_substrings requires exactly two strings")

        s1, s2 = sequences
        common = set()
        s2_length = len(s2)

        for start in range(s2_length):
            for end in range(start + 2, s2_length + 1):  # substrings of length >= 2
                substr = s2[start:end]
                if substr not in s1:
                    break
                common.add(substr)

        return list(common)

    def __call__(self, sequences: list[str]) -> list[str]:
        if (
            not isinstance(sequences, list) and not isinstance(sequences, tuple)
        ) or not all(isinstance(s, str) for s in sequences):
            raise TypeError("longest_common_substring expects a list of strings")
        if len(sequences) < 2:
            raise ValueError("Provide at least two sequences")

        # Generate substrings from first and last strings
        motifs = self._common_substrings([sequences[0], sequences[-1]])
        if not motifs:
            return [""]

        shared = [motif for motif in motifs if all(motif in seq for seq in sequences)]
        if not shared:
            return [""]

        max_len = len(max(shared, key=len))
        return [m for m in shared if len(m) == max_len]

    def align(self, sequences: list[str]) -> list[str]:
        return self(sequences)

    def distance(self, sequences: list[str]) -> int:
        return len(max(sequences, key=len)) - self.similarity(sequences)

    def similarity(self, sequences: list[str]) -> int:
        lcsub = self(sequences)
        if not lcsub:
            return 0
        return len(lcsub[0])

    def normalized_distance(self, sequences: list[str]) -> float:
        return 1 - self.normalized_similarity(sequences)

    def normalized_similarity(self, sequences: list[str]) -> float:
        lcsub = self(sequences)
        if not lcsub:
            return 0
        return len(lcsub[0]) / len(max(sequences, key=len))


feng_doolittle = FengDoolittle()
longest_common_substring_msa = LongestCommonSubstringMSA()


if __name__ == "__main__":
    main()
