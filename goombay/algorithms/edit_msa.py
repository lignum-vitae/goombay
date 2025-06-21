try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
except ImportError:
    raise ImportError("Please pip install all dependencies from requirements.txt!")

# internal dependencies
from goombay.algorithms.edit import (
    NeedlemanWunsch,
    LowranceWagner,
    WagnerFischer,
    WatermanSmithBeyer,
    Gotoh,
    Hirschberg,
    Jaro,
    JaroWinkler,
)


def main():
    seq1 = "HOUSEOFCARDSFALLDOWN"
    seq2 = "HOUSECARDFALLDOWN"
    seq3 = "FALLDOWN"


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


feng_doolittle = FengDoolittle()

if __name__ == "__main__":
    main()
