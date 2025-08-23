# standard library
from abc import ABC, abstractmethod

# external dependencies
from numpy import float64
from numpy.typing import NDArray


class GlobalBase(ABC):
    @abstractmethod
    def __call__(
        self, query_seq: str, subject_seq: str
    ) -> tuple[NDArray[float64], NDArray[float64]]:
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
        if query_seq == subject_seq:
            return 1.0
        if not query_seq or not subject_seq:
            return 0.0

        raw_score = self.similarity(query_seq, subject_seq)
        if self.has_sub_mat:
            max_possible = 0
            min_possible = 0
            for q, s in zip(query_seq, subject_seq):
                q_match = self.sub_mat[q][q]
                s_match = self.sub_mat[s][s]
                qs_match = self.sub_mat[q][s]
                candidates = (q_match, s_match, qs_match)
                max_possible += max(candidates)
                min_possible += min(candidates)
        else:
            max_len = len(max(query_seq, subject_seq, key=len))
            max_possible = max_len * self.match
            min_possible = -max_len * self.mismatch
        score_range = max_possible + abs(min_possible)
        return (raw_score + abs(min_possible)) / score_range

    def align(
        self, query_seq: str, subject_seq: str, all_alignments: bool = False
    ) -> str | list[str]:
        _, pointer_matrix = self(query_seq, subject_seq)

        qs = [x.upper() for x in query_seq]
        ss = [x.upper() for x in subject_seq]
        i, j = len(qs), len(ss)
        aligned = []
        stack = [([""], [""], i, j)]

        # looks for match/mismatch/gap starting from bottom right of matrix
        while stack:
            qs_align, ss_align, i, j = stack.pop()
            if i <= 0 and j <= 0:
                qs_aligned = "".join(qs_align[::-1])
                ss_aligned = "".join(ss_align[::-1])
                aligned.append(f"{qs_aligned}\n{ss_aligned}")
                continue
            if pointer_matrix[i, j] in [2, 2 + 3, 2 + 4, 2 + 3 + 4]:
                # appends match/mismatch then moves to the cell diagonally up and to the left
                stack.append(
                    (qs_align + [qs[i - 1]], ss_align + [ss[j - 1]], i - 1, j - 1)
                )
                if not all_alignments:
                    continue
            if pointer_matrix[i, j] in [3, 3 + 2, 3 + 4, 3 + 2 + 4]:
                # appends gap and accompanying nucleotide, then moves to the cell above
                stack.append((qs_align + [qs[i - 1]], ss_align + ["-"], i - 1, j))
                if not all_alignments:
                    continue
            if pointer_matrix[i, j] in [4, 4 + 2, 4 + 3, 4 + 2 + 3]:
                # appends gap and accompanying nucleotide, then moves to the cell to the left
                stack.append((qs_align + ["-"], ss_align + [ss[j - 1]], i, j - 1))
                if not all_alignments:
                    continue

        if not all_alignments:
            return aligned[0]
        return aligned


class LocalBase(ABC):
    @abstractmethod
    def __call__(self, query_seq: str, subject_seq: str) -> NDArray[float64]:
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
