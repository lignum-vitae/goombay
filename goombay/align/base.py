# standard library
from abc import ABC, abstractmethod

# external dependencies
from numpy.typing import NDArray


class GlobalBase(ABC):
    @abstractmethod
    def __call__(self, query_seq: str, subject_seq: str):
        pass

    def matrix(self, query_seq: str, subject_seq: str) -> list[list[float]]:
        matrix, _ = self(query_seq, subject_seq)
        return matrix

    def distance(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq and not subject_seq:
            return 0.0
        if not query_seq or not subject_seq:
            return float(len(query_seq or subject_seq)) * self.gap

        raw_sim = self.similarity(query_seq, subject_seq)
        max_possible = max(len(query_seq), len(subject_seq)) * self.match
        return max_possible - abs(raw_sim)

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq and not subject_seq:
            return 1.0
        matrix, _ = self(query_seq, subject_seq)
        return matrix[matrix.shape[0] - 1, matrix.shape[1] - 1]

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return 1 - self.normalized_similarity(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        raw_score = self.similarity(query_seq, subject_seq)
        max_len = len(max(query_seq, subject_seq, key=len))
        max_possible = max_len * self.match
        min_possible = -max_len * self.mismatch
        score_range = max_possible - min_possible
        return (raw_score - min_possible) / score_range

    def align(self, query_seq: str, subject_seq: str) -> str:
        _, pointer_matrix = self(query_seq, subject_seq)

        qs = [x.upper() for x in query_seq]
        ss = [x.upper() for x in subject_seq]
        i, j = len(qs), len(ss)
        qs_align, ss_align = [], []

        # looks for match/mismatch/gap starting from bottom right of matrix
        while i > 0 or j > 0:
            if pointer_matrix[i, j] in [2, 5, 6, 9]:
                # appends match/mismatch then moves to the cell diagonally up and to the left
                qs_align.append(qs[i - 1])
                ss_align.append(ss[j - 1])
                i -= 1
                j -= 1
            elif pointer_matrix[i, j] in [3, 5, 7, 9]:
                # appends gap and accompanying nucleotide, then moves to the cell above
                ss_align.append("-")
                qs_align.append(qs[i - 1])
                i -= 1
            elif pointer_matrix[i, j] in [4, 6, 7, 9]:
                # appends gap and accompanying nucleotide, then moves to the cell to the left
                ss_align.append(ss[j - 1])
                qs_align.append("-")
                j -= 1

        qs_align = "".join(qs_align[::-1])
        ss_align = "".join(ss_align[::-1])

        return f"{qs_align}\n{ss_align}"


class LocalBase(ABC):
    @abstractmethod
    def __call__(self, query_seq: str, subject_seq: str):
        pass

    def matrix(self, query_seq: str, subject_seq: str) -> NDArray:
        """Return alignment matrix"""
        return self(query_seq, subject_seq)

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        """Calculate similarity score"""
        if not query_seq and not subject_seq:
            return 1.0
        if not query_seq or not subject_seq:
            return 0.0
        if len(query_seq) == 1 and len(subject_seq) == 1 and query_seq == subject_seq:
            return 1.0
        matrix = self(query_seq, subject_seq)
        return matrix.max() if matrix.max() > 1 else 0.0

    def distance(self, query_seq: str, subject_seq: str) -> float:
        query_length = len(query_seq)
        subject_length = len(subject_seq)
        if not query_seq and not subject_seq:
            return 0.0
        if not query_seq or not subject_seq:
            return max(query_length, subject_length)

        matrix = self(query_seq, subject_seq)
        sim_AB = matrix.max()
        max_score = self.match * max(query_length, subject_length)
        return max_score - sim_AB

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        """Calculate normalized similarity between 0 and 1"""
        if not query_seq and not subject_seq:
            return 1.0
        if not query_seq or not subject_seq:
            return 0.0
        if len(query_seq) == 1 and len(subject_seq) == 1 and query_seq == subject_seq:
            return 1.0
        matrix = self(query_seq, subject_seq)
        best_score = matrix.max()
        return best_score / min(len(query_seq), len(subject_seq))

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        """Calculate normalized distance between 0 and 1"""
        if not query_seq and not subject_seq:
            return 0.0
        if not query_seq or not subject_seq:
            return 1.0
        return 1.0 - self.normalized_similarity(query_seq, subject_seq)
