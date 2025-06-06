# internal dependencies
from goombay.algorithms.editdistance import (
    Needleman_Wunsch,
    Lowrance_Wagner,
    Wagner_Fischer,
    Waterman_Smith_Beyer,
    Gotoh,
    Hirschberg,
    Jaro,
    Jaro_Winkler,
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

    # print(feng_doolittle.align(seq1, seq2, seq3))
    print(feng_doolittle.pairwise)
    print(Feng_Doolittle.supported_pairwise_algs())


class Feng_Doolittle:
    supported_pairwise = {
        "needleman_wunsch": Needleman_Wunsch,
        "jaro": Jaro,
        "jaro_winkler": Jaro_Winkler,
        "gotoh": Gotoh,
        "wagner_fischer": Wagner_Fischer,
        "waterman_smith_beyer": Waterman_Smith_Beyer,
        "hirschberg": Hirschberg,
        "lowrance_wagner": Lowrance_Wagner,
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


feng_doolittle = Feng_Doolittle()

if __name__ == "__main__":
    main()
