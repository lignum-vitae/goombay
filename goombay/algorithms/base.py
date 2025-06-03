import math
from numpy.typing import NDArray

class GLOBALBASE:
    def matrix(self, querySequence: str, subjectSequence: str) -> list[list[float]]:
        matrix, _ = self(querySequence, subjectSequence)
        return matrix

    def distance(self, querySequence: str, subjectSequence: str) -> float:
        if not querySequence and not subjectSequence:
            return 0.0
        if not querySequence or not subjectSequence:
            return float(len(querySequence or subjectSequence)) * self.gap_penalty

        raw_sim = self.similarity(querySequence, subjectSequence)
        max_possible = max(len(querySequence), len(subjectSequence)) * self.match_score
        return max_possible - abs(raw_sim)

    def similarity(self, querySequence: str, subjectSequence: str) -> float:
        if not querySequence and not subjectSequence:
            return 1.0
        matrix, _ = self(querySequence, subjectSequence)
        return matrix[matrix.shape[0]-1,matrix.shape[1]-1]

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        return 1 -  self.normalized_similarity(querySequence, subjectSequence)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        raw_score = self.similarity(querySequence, subjectSequence)
        max_len = len(max(querySequence, subjectSequence, key=len))
        max_possible = max_len * self.match_score
        min_possible = -max_len * self.mismatch_penalty
        score_range = max_possible - min_possible
        return (raw_score - min_possible) / score_range


    def align(self, querySequence: str, subjectSequence: str) -> str:
        _, pointerMatrix = self(querySequence, subjectSequence)

        qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        i, j = len(qs), len(ss)
        queryAlign, subjectAlign = [], []

        while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
          if pointerMatrix[i,j] in [2, 5, 6, 9]:
              #appends match/mismatch then moves to the cell diagonally up and to the left
              queryAlign.append(qs[i-1])
              subjectAlign.append(ss[j-1])
              i -= 1
              j -= 1
          elif pointerMatrix[i,j] in [3, 5, 7, 9]:
              #appends gap and accompanying nucleotide, then moves to the cell above
              subjectAlign.append('-')
              queryAlign.append(qs[i-1])
              i -= 1
          elif pointerMatrix[i,j] in [4, 6, 7, 9]:
              #appends gap and accompanying nucleotide, then moves to the cell to the left
              subjectAlign.append(ss[j-1])
              queryAlign.append('-')
              j -= 1

        queryAlign = "".join(queryAlign[::-1])
        subjectAlign = "".join(subjectAlign[::-1])

        return f"{queryAlign}\n{subjectAlign}"

class LOCALBASE:
    def matrix(self, querySequence: str, subjectSequence: str) -> NDArray:
        """Return alignment matrix"""
        matrix = self(querySequence, subjectSequence)
        return matrix

    def similarity(self, querySequence: str, subjectSequence: str) -> float:
        """Calculate similarity score"""
        matrix = self(querySequence, subjectSequence)
        return matrix.max()

    def _compute_three_way_similarity(self, querySequence: str, subjectSequence: str) -> tuple[float, float, float]:
        """Compute similarity between sequences and their self-similarities efficiently."""
        if not querySequence and not subjectSequence:
            return 0.0, 0.0, 0.0
        if not querySequence or not subjectSequence:
            return 0.0, 0.0, 0.0

        if querySequence == subjectSequence:
            matrix = self(querySequence, querySequence)
            if isinstance(matrix, tuple):
                matrix = matrix[0]
            sim = matrix.max()
            return sim, sim, sim

        return self._compute_all_matrices(querySequence, subjectSequence)

    def _compute_all_matrices(self, querySequence: str, subjectSequence: str) -> tuple[float, float, float]:
        """Template method for computing all three alignments."""
        # Initialize matrices
        matrices_A = self._init_alignment_matrices(len(querySequence), len(querySequence))
        matrices_B = self._init_alignment_matrices(len(subjectSequence), len(subjectSequence))
        matrices_AB = self._init_alignment_matrices(len(querySequence), len(subjectSequence))

        # Fill matrices
        for i in range(1, len(querySequence) + 1):
            # Fill A matrices (querySequence self-alignment)
            for j in range(1, len(querySequence) + 1):
                self._fill_cell(matrices_A, i, j, querySequence[i-1], querySequence[j-1])

            # Fill B matrices (subjectSequence self-alignment)
            if i <= len(subjectSequence):
                for j in range(1, len(subjectSequence) + 1):
                    self._fill_cell(matrices_B, i, j, subjectSequence[i-1], subjectSequence[j-1])

            # Fill AB matrices (querySequence-subjectSequence alignment)
            for j in range(1, len(subjectSequence) + 1):
                self._fill_cell(matrices_AB, i, j, querySequence[i-1], subjectSequence[j-1])

        return self._get_max_scores(matrices_A, matrices_B, matrices_AB)

    def distance(self, querySequence: str, subjectSequence: str) -> float:
        """Calculate a proper metric distance based on local alignment score.

        Uses the formula: d(x,y) = -ln(sim_AB / sqrt(sim_A * sim_B))
        This ensures the triangle inequality property.
        """
        if not querySequence and not subjectSequence:
            return 0.0
        if not querySequence or not subjectSequence:
            return max(len(querySequence), len(subjectSequence))

        sim_A, sim_B, sim_AB = self._compute_three_way_similarity(querySequence, subjectSequence)
        if sim_AB == 0:
            return max(len(querySequence), len(subjectSequence))
        return -math.log(sim_AB / math.sqrt(sim_A * sim_B))

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        """Calculate normalized similarity between 0 and 1"""
        if not querySequence and not subjectSequence:
            return 1.0
        if not querySequence or not subjectSequence:
            return 0.0
        similarity = self.similarity(querySequence, subjectSequence)
        opt_score = min(len(querySequence), len(subjectSequence)) * self.match_score
        return similarity/opt_score

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        """Calculate normalized distance between 0 and 1"""
        if not querySequence and not subjectSequence:
            return 0.0
        if not querySequence or not subjectSequence:
            return 1.0
        return 1.0 - self.normalized_similarity(querySequence, subjectSequence)
