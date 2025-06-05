import math
from numpy.typing import NDArray

class GlobalBase:
    def matrix(self, query_sequence: str, subject_sequence: str) -> list[list[float]]:
        matrix, _ = self(query_sequence, subject_sequence)
        return matrix

    def distance(self, query_sequence: str, subject_sequence: str) -> float:
        if not query_sequence and not subject_sequence:
            return 0.0
        if not query_sequence or not subject_sequence:
            return float(len(query_sequence or subject_sequence)) * self.gap_penalty

        raw_sim = self.similarity(query_sequence, subject_sequence)
        max_possible = max(len(query_sequence), len(subject_sequence)) * self.match_score
        return max_possible - abs(raw_sim)

    def similarity(self, query_sequence: str, subject_sequence: str) -> float:
        if not query_sequence and not subject_sequence:
            return 1.0
        matrix, _ = self(query_sequence, subject_sequence)
        return matrix[matrix.shape[0]-1,matrix.shape[1]-1]

    def normalized_distance(self, query_sequence: str, subject_sequence: str) -> float:
        return 1 -  self.normalized_similarity(query_sequence, subject_sequence)

    def normalized_similarity(self, query_sequence: str, subject_sequence: str) -> float:
        raw_score = self.similarity(query_sequence, subject_sequence)
        max_len = len(max(query_sequence, subject_sequence, key=len))
        max_possible = max_len * self.match_score
        min_possible = -max_len * self.mismatch_penalty
        score_range = max_possible - min_possible
        return (raw_score - min_possible) / score_range


    def align(self, query_sequence: str, subject_sequence: str) -> str:
        _, pointer_matrix = self(query_sequence, subject_sequence)

        qs, ss = [x.upper() for x in query_sequence], [x.upper() for x in subject_sequence]
        i, j = len(qs), len(ss)
        query_align, subject_align = [], []

        while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
          if pointer_matrix[i,j] in [2, 5, 6, 9]:
              #appends match/mismatch then moves to the cell diagonally up and to the left
              query_align.append(qs[i-1])
              subject_align.append(ss[j-1])
              i -= 1
              j -= 1
          elif pointer_matrix[i,j] in [3, 5, 7, 9]:
              #appends gap and accompanying nucleotide, then moves to the cell above
              subject_align.append('-')
              query_align.append(qs[i-1])
              i -= 1
          elif pointer_matrix[i,j] in [4, 6, 7, 9]:
              #appends gap and accompanying nucleotide, then moves to the cell to the left
              subject_align.append(ss[j-1])
              query_align.append('-')
              j -= 1

        query_align = "".join(query_align[::-1])
        subject_align = "".join(subject_align[::-1])

        return f"{query_align}\n{subject_align}"

class LocalBase:
    def matrix(self, query_sequence: str, subject_sequence: str) -> NDArray:
        """Return alignment matrix"""
        matrix = self(query_sequence, subject_sequence)
        return matrix

    def similarity(self, query_sequence: str, subject_sequence: str) -> float:
        """Calculate similarity score"""
        matrix = self(query_sequence, subject_sequence)
        return matrix.max()

    def _compute_three_way_similarity(self, query_sequence: str, subject_sequence: str) -> tuple[float, float, float]:
        """Compute similarity between sequences and their self-similarities efficiently."""
        if not query_sequence or not subject_sequence:
            return 0.0, 0.0, 0.0

        if query_sequence == subject_sequence:
            matrix = self(query_sequence, query_sequence)
            if isinstance(matrix, tuple):
                matrix = matrix[0]
            sim = matrix.max()
            return sim, sim, sim

        return self._compute_all_matrices(query_sequence, subject_sequence)

    def _compute_all_matrices(self, query_sequence: str, subject_sequence: str) -> tuple[float, float, float]:
        """Template method for computing all three alignments."""
        query_length = len(query_sequence)
        subject_length = len(subject_sequence)

        # Initialize matrices
        matrices_A = self._init_alignment_matrices(len(query_sequence), len(query_sequence))
        matrices_B = self._init_alignment_matrices(len(subject_sequence), len(subject_sequence))
        matrices_AB = self._init_alignment_matrices(len(query_sequence), len(subject_sequence))

        # Fill matrices
        for i in range(1, query_length + 1):
            # Fill A matrices (query_sequence self-alignment)
            for j in range(1, query_length + 1):
                self._fill_cell(matrices_A, i, j, query_sequence[i-1], query_sequence[j-1])

            # Fill B matrices (subject_sequence self-alignment)
            if i <= subject_length:
                for j in range(1, subject_length + 1):
                    self._fill_cell(matrices_B, i, j, subject_sequence[i-1], subject_sequence[j-1])

            # Fill AB matrices (query_sequence-subject_sequence alignment)
            for j in range(1, subject_length + 1):
                self._fill_cell(matrices_AB, i, j, query_sequence[i-1], subject_sequence[j-1])

        return self._get_max_scores(matrices_A, matrices_B, matrices_AB)

    def distance(self, query_sequence: str, subject_sequence: str) -> float:
        """Calculate a proper metric distance based on local alignment score.

        Uses the formula: d(x,y) = -ln(sim_AB / sqrt(sim_A * sim_B))
        This ensures the triangle inequality property.
        """
        query_length = len(query_sequence)
        subject_length = len(subject_sequence)

        if not query_sequence and not subject_sequence:
            return 0.0
        if not query_sequence or not subject_sequence:
            return max(query_length, subject_length)

        sim_A, sim_B, sim_AB = self._compute_three_way_similarity(query_sequence, subject_sequence)
        if sim_AB == 0:
            return max(query_length, subject_length)
        return -math.log(sim_AB / math.sqrt(sim_A * sim_B))

    def normalized_similarity(self, query_sequence: str, subject_sequence: str) -> float:
        """Calculate normalized similarity between 0 and 1"""
        if not query_sequence and not subject_sequence:
            return 1.0
        if not query_sequence or not subject_sequence:
            return 0.0
        similarity = self.similarity(query_sequence, subject_sequence)
        opt_score = min(len(query_sequence), len(subject_sequence)) * self.match_score
        return similarity/opt_score

    def normalized_distance(self, query_sequence: str, subject_sequence: str) -> float:
        """Calculate normalized distance between 0 and 1"""
        if not query_sequence and not subject_sequence:
            return 0.0
        if not query_sequence or not subject_sequence:
            return 1.0
        return 1.0 - self.normalized_similarity(query_sequence, subject_sequence)
