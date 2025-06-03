#internal dependencies
from goombay.algorithms.base import GLOBALBASE as _GLOBALBASE, LOCALBASE as _LOCALBASE

try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
except ImportError:
    raise ImportError("Please pip install all dependencies from requirements.txt!")

def main():
    """
    qqs = "HOLYWATERISABLESSING"
    sss = ["HOLYWATERBLESSING", "HOLYERISSING", "HOLYWATISSSI", "HWATISBLESSING", "HOLYWATISSS"]

    for i in range(len(sss)):
        print(waterman_smith_beyer.align(qqs, sss[i]))
        print()
    print(waterman_smith_beyer.matrix("TRATE", "TRACE"))
    """
    query = "AGCT"
    subject = "TACG"
    print(needleman_wunsch.distance(query, subject))
    print(needleman_wunsch.similarity(query, subject))
    print(needleman_wunsch.normalized_similarity(query, subject))
    print(needleman_wunsch.normalized_distance(query, subject))
    print(needleman_wunsch.matrix(query, subject))
    print(needleman_wunsch.align(query, subject))

class Wagner_Fischer(_GLOBALBASE): #Levenshtein Distance
    def __init__(self)->None:
        self.gap_penalty = 1
        self.substitution_cost = 1

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64]]:
        qs,ss = [""], [""]
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [n for n in range(len(qs))]
        self.alignment_score[0,:] = [n for n in range(len(ss))]

        for i in range(1, len(qs)):
          for j in range(1, len(ss)):
              substitution_cost = 0
              if qs[i] != ss[j]:
                  substitution_cost = self.substitution_cost
              match = self.alignment_score[i-1][j-1] + substitution_cost
              ugap = self.alignment_score[i-1][j] + self.gap_penalty
              lgap = self.alignment_score[i][j-1] + self.gap_penalty

              tmin = min(match, lgap, ugap)

              self.alignment_score[i][j] = tmin #lowest value is best choice
              if match == tmin: #matrix for traceback based on results from scoring matrix
                  self.pointer[i,j] += 2
              if ugap == tmin:
                  self.pointer[i,j] += 3
              if lgap == tmin:
                  self.pointer[i,j] += 4
        return self.alignment_score, self.pointer

    def distance(self, querySequence: str, subjectSequence: str)->float:
        matrix, _ = self(querySequence, subjectSequence)
        return float(matrix[-1, -1])

    def similarity(self, querySequence: str, subjectSequence: str) -> float:
        if not querySequence and not subjectSequence:
            return 1.0
        sim = max(len(querySequence), len(subjectSequence)) - self.distance(querySequence, subjectSequence)
        return max(0, sim)

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        if not querySequence and not subjectSequence:
            return 0.0
        if not querySequence or not subjectSequence:
            return 1.0
        max_len = max(len(str(querySequence)), len(str(subjectSequence)))
        max_dist = max_len
        return (self.distance(querySequence, subjectSequence) / max_dist)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        return 1.0 - self.normalized_distance(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
        return super().matrix(querySequence, subjectSequence)

    def align(self, querySequence: str, subjectSequence: str)->str:
        _, pointerMatrix = self(querySequence, subjectSequence)

        qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        i, j = len(qs), len(ss)
        queryAlign, subjectAlign = [], []
        while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
          if pointerMatrix[i,j] in [2, 5, 6, 10, 9, 13, 14, 17]:
              #appends match/mismatch then moves to the cell diagonally up and to the left
              queryAlign.append(qs[i-1])
              subjectAlign.append(ss[j-1])
              i -= 1
              j -= 1
          elif pointerMatrix[i,j] in [3, 5, 7, 11, 9, 13, 15, 17]:
              #appends gap and accompanying nucleotide, then moves to the cell above
              subjectAlign.append('-')
              queryAlign.append(qs[i-1])
              i -= 1
          elif pointerMatrix[i,j] in [4, 6, 7, 12, 9, 14, 15, 17]:
              #appends gap and accompanying nucleotide, then moves to the cell to the left
              subjectAlign.append(ss[j-1])
              queryAlign.append('-')
              j -= 1

        queryAlign = "".join(queryAlign[::-1])
        subjectAlign = "".join(subjectAlign[::-1])
        return f"{queryAlign}\n{subjectAlign}"

class Lowrance_Wagner(_GLOBALBASE): #Damerau-Levenshtein distance
    def __init__(self)->None:
        self.gap_penalty = 1
        self.substitution_cost = 1
        self.transposition_cost = 1

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64]]:
        qs,ss = [""], [""]
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [n for n in range(len(qs))]
        self.alignment_score[0,:] = [n for n in range(len(ss))]

        for i in range(1, len(qs)):
            for j in range(1, len(ss)):
                substitution_cost = 0
                if qs[i] != ss[j]:
                    substitution_cost = self.substitution_cost
                match = self.alignment_score[i-1][j-1] + substitution_cost
                ugap = self.alignment_score[i-1][j] + self.gap_penalty
                lgap = self.alignment_score[i][j-1] + self.gap_penalty
                trans = self.alignment_score[i-2][j-2] + 1 if qs[i] == ss[j-1] and ss[j] == qs[i-1] else float('inf')
                tmin = min(match, lgap, ugap, trans)

                self.alignment_score[i][j] = tmin #lowest value is best choice
                if match == tmin: #matrix for traceback based on results from scoring matrix
                    self.pointer[i,j] += 2
                if ugap == tmin:
                    self.pointer[i,j] += 3
                if lgap == tmin:
                    self.pointer[i,j] += 4
                if trans == tmin:
                    self.pointer[i,j] += 8
        return self.alignment_score, self.pointer

    def distance(self, querySequence: str, subjectSequence: str)->float:
        matrix, _ = self(querySequence, subjectSequence)
        return float(matrix[-1, -1])

    def similarity(self, querySequence: str, subjectSequence: str) -> float:
        if not querySequence and not subjectSequence:
            return 1.0
        sim = max(len(querySequence), len(subjectSequence)) - self.distance(querySequence, subjectSequence)
        return max(0, sim)

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        if not querySequence and not subjectSequence:
            return 0.0
        if not querySequence or not subjectSequence:
            return 1.0
        max_len = max(len(str(querySequence)), len(str(subjectSequence)))
        max_dist = max_len
        return (self.distance(querySequence, subjectSequence) / max_dist)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        return 1.0 - self.normalized_distance(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
        return super().matrix(querySequence, subjectSequence)

    def align(self, querySequence: str, subjectSequence: str)->str:
        if not querySequence and not subjectSequence:
            return "\n"
        if not querySequence:
            return f"{'-' * len(subjectSequence)}\n{subjectSequence}"
        if not subjectSequence:
            return f"{querySequence}\n{'-' * len(querySequence)}"

        _, pointerMatrix = self(querySequence, subjectSequence)

        qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        i, j = len(qs), len(ss)
        queryAlign, subjectAlign = [], []
        while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
          if pointerMatrix[i,j] in [2, 5, 6, 10, 9, 13, 14, 17]:
              #appends match/mismatch then moves to the cell diagonally up and to the left
              queryAlign.append(qs[i-1])
              subjectAlign.append(ss[j-1])
              i -= 1
              j -= 1
          elif pointerMatrix[i,j] in [8, 10, 11, 12, 13, 14, 15, 17]:
              queryAlign.extend([qs[i-1],qs[i-2]])
              subjectAlign.extend([ss[j-1],ss[j-2]])
              i -= 2
              j-= 2
          elif pointerMatrix[i,j] in [3, 5, 7, 11, 9, 13, 15, 17]:
              #appends gap and accompanying nucleotide, then moves to the cell above
              subjectAlign.append('-')
              queryAlign.append(qs[i-1])
              i -= 1
          elif pointerMatrix[i,j] in [4, 6, 7, 12, 9, 14, 15, 17]:
              #appends gap and accompanying nucleotide, then moves to the cell to the left
              subjectAlign.append(ss[j-1])
              queryAlign.append('-')
              j -= 1

        queryAlign = "".join(queryAlign[::-1])
        subjectAlign = "".join(subjectAlign[::-1])
        return f"{queryAlign}\n{subjectAlign}"

class Hamming():
    def _check_inputs(self, querySequence: str|int, subjectSequence: str|int) -> None:
        if not isinstance(querySequence, (str, int)) or not isinstance(subjectSequence, (str, int)):
            raise TypeError("Sequences must be strings or integers")
        if type(querySequence) != type(subjectSequence):
            raise TypeError("Sequences must be of the same type (both strings or both integers)")
        if len(str(querySequence)) != len(str(subjectSequence)) and not isinstance(querySequence, int):
            raise ValueError("Sequences must be of equal length")

    def __call__(self, querySequence: str|int, subjectSequence: str|int)->tuple[int,list[int]]:
        self._check_inputs(querySequence, subjectSequence)
        if isinstance(querySequence, int) and isinstance(subjectSequence, int):
            qs, ss = bin(querySequence)[2:], bin(subjectSequence)[2:]
            # Pad with leading zeros to make equal length
            max_len = max(len(qs), len(ss))
            qs = qs.zfill(max_len)
            ss = ss.zfill(max_len)
        else:
            qs,ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]

        if len(qs) == 1 and len(ss) == 1:
            dist = 1 if qs != ss else 0
            dist_array = [dist]
            return dist, dist_array

        dist = 0
        dist_array = []
        for i, char in enumerate(qs):
            if char != ss[i]:
                dist += 1
                dist_array.append(1)
                continue
            dist_array.append(0)

        dist += len(ss)-len(qs)
        dist_array.extend([1]*(len(ss)-len(qs)))
        return dist, dist_array

    def distance(self, querySequence: str|int, subjectSequence: str|int)->int:
        self._check_inputs(querySequence, subjectSequence)
        if isinstance(querySequence, int) and isinstance(subjectSequence, int):
            qs, ss = int(querySequence), int(subjectSequence)
            return bin(qs ^ ss).count("1")
        if len(querySequence) == len(subjectSequence) == 0:
            return 0
        qs,ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        query = set([(x, y) for (x, y) in enumerate(qs)])
        subject = set([(x, y) for (x, y) in enumerate(ss)])
        qs,sq = query-subject, subject-query
        dist = max(map(len,[qs,sq]))
        return dist 

    def similarity(self, querySequence: str|int, subjectSequence: str|int)->int:
        self._check_inputs(querySequence, subjectSequence)
        if isinstance(querySequence, int) and isinstance(subjectSequence, int):
            qs, ss = int(querySequence), int(subjectSequence)
            return bin(qs & ss).count("1")
        if len(querySequence) == len(subjectSequence) == 0:
            return 1
        qs,ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        query = set([(x, y) for (x, y) in enumerate(qs)])
        subject = set([(x, y) for (x, y) in enumerate(ss)])
        qs,sq = query-subject, subject-query
        sim = max(map(len, [querySequence, subjectSequence])) - max(map(len,[qs,sq]))
        return sim

    def normalized_distance(self, querySequence, subjectSequence) -> float:
        return self.distance(querySequence, subjectSequence)/len(querySequence)

    def normalized_similarity(self, querySequence, subjectSequence) -> float:
        return 1 - self.normalized_distance(querySequence, subjectSequence)

    def binary_distance_array(self, querySequence: str, subjectSequence: str)->list[int]:
        self._check_inputs(querySequence, subjectSequence)
        _, distarray = self(querySequence, subjectSequence)
        return distarray

    def binary_similarity_array(self, querySequence: str, subjectSequence: str)->list[int]:
        self._check_inputs(querySequence, subjectSequence)
        _, distarray = self(querySequence, subjectSequence)
        simarray = [1 if num == 0 else 0 for num in distarray]
        return simarray

    def matrix(self, qs: str, ss: str) -> None:
        return None

    def align(self, querySequence: str|int, subjectSequence: str|int)->str:
        self._check_inputs(querySequence, subjectSequence)
        if isinstance(querySequence, int) and isinstance(subjectSequence, int):
            qs, ss = int(querySequence), int(subjectSequence)
            return f"{bin(qs)}\n{bin(ss)}"
        return f"{querySequence}\n{subjectSequence}"

class Needleman_Wunsch(_GLOBALBASE):
    def __init__(self, match_score:int = 2, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64]]:
        qs,ss = [""], [""] 
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [-n*self.gap_penalty for n in range(len(qs))]
        self.alignment_score[0,:] = [-n*self.gap_penalty for n in range(len(ss))]

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
              if i == 0 or j == 0:
                  #keeps first row and column consistent throughout all calculations
                  continue
              if query_char == subject_char:
                  match = self.alignment_score[i-1][j-1] + self.match_score
              else:
                  match = self.alignment_score[i-1][j-1] - self.mismatch_penalty
              ugap = self.alignment_score[i-1][j] - self.gap_penalty
              lgap = self.alignment_score[i][j-1] - self.gap_penalty
              tmax = max(match, lgap, ugap)

              self.alignment_score[i][j] = tmax #highest value is best choice
              if match == tmax: #matrix for traceback based on results from scoring matrix
                  self.pointer[i,j] += 2
              if ugap == tmax:
                  self.pointer[i,j] += 3
              if lgap == tmax:
                  self.pointer[i,j] += 4
        return self.alignment_score, self.pointer

    def distance(self, querySequence: str, subjectSequence: str)->float:
        return super().distance(querySequence, subjectSequence)

    def similarity(self, querySequence: str, subjectSequence: str)->float:
        return super().similarity(querySequence, subjectSequence)

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_distance(querySequence, subjectSequence)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_similarity(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
        return super().matrix(querySequence, subjectSequence)

    def align(self, querySequence: str, subjectSequence: str) -> str:
        return super().align(querySequence, subjectSequence)

class Waterman_Smith_Beyer(_GLOBALBASE):
    def __init__(self, match_score:int = 2, mismatch_penalty:int = 1, new_gap_penalty:int = 4, continue_gap_penalty:int = 1)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.new_gap_penalty = new_gap_penalty
        self.continue_gap_penalty = continue_gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64], NDArray[float64]]:
        qs,ss = [""], [""] 
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [-self.new_gap_penalty + -n * self.continue_gap_penalty for n in range(len(qs))]
        self.alignment_score[0,:] = [-self.new_gap_penalty + -n * self.continue_gap_penalty for n in range(len(ss))]
        self.alignment_score[0][0] = 0

        for i, subject in enumerate(qs):
            for j, query in enumerate(ss):
                if i == 0 or j == 0:
                    #keeps first row and column consistent throughout all calculations
                    continue
                if subject == query: 
                    matchScore = self.alignment_score[i-1][j-1] + self.match_score
                else:
                    matchScore = self.alignment_score[i-1][j-1] - self.mismatch_penalty
                #both gaps defaulted to continue gap penalty
                ugapScore = self.alignment_score[i-1][j] - self.continue_gap_penalty
                lgapScore = self.alignment_score[i][j-1] - self.continue_gap_penalty
                #if cell before i-1 or j-1 is gap, then this is a gap continuation
                if self.alignment_score[i-1][j] != (self.alignment_score[i-2][j]) - self.new_gap_penalty - self.continue_gap_penalty:
                    ugapScore -= self.new_gap_penalty
                if self.alignment_score[i][j-1] != (self.alignment_score[i][j-2]) - self.new_gap_penalty - self.continue_gap_penalty:
                    lgapScore -= self.new_gap_penalty
                tmax = max(matchScore, lgapScore, ugapScore)

                self.alignment_score[i][j] = tmax #highest value is best choice
                #matrix for traceback based on results from scoring matrix
                if matchScore == tmax: 
                    self.pointer[i,j] += 2
                elif ugapScore == tmax:
                    self.pointer[i,j] += 3
                elif lgapScore == tmax:
                    self.pointer[i,j] += 4
        return self.alignment_score, self.pointer

    def distance(self, querySequence: str, subjectSequence: str)->float:
        return super().distance(querySequence, subjectSequence)

    def similarity(self, querySequence: str, subjectSequence: str)->float:
        return super().similarity(querySequence, subjectSequence)

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_distance(querySequence, subjectSequence)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_similarity(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
        return super().matrix(querySequence, subjectSequence)

    def align(self, querySequence: str, subjectSequence: str) -> str:
        return super().align(querySequence, subjectSequence)

class Gotoh(_GLOBALBASE):
    def __init__(self, match_score:int = 2, mismatch_penalty:int = 1, new_gap_penalty:int = 2, continue_gap_penalty: int = 1)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.new_gap_penalty = new_gap_penalty
        self.continue_gap_penalty = continue_gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64],NDArray[float64],NDArray[float64]]:
        qs,ss = [""], [""] 
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.D = numpy.full((len(qs),len(ss)), -numpy.inf)
        self.P = numpy.full((len(qs), len(ss)), -numpy.inf)
        self.P[:,0] = 0
        self.Q = numpy.full((len(qs), len(ss)), -numpy.inf)
        self.Q[0,:] = 0
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.D[0, 0] = 0
        # Initialize first column (vertical gaps)
        for i in range(1, len(qs)):
            self.D[i, 0] = -(self.new_gap_penalty + (i) * self.continue_gap_penalty)
        # Initialize first row (horizontal gaps)
        for j in range(1, len(ss)):
            self.D[0, j] = -(self.new_gap_penalty + (j) * self.continue_gap_penalty)

        for i in range(1, len(qs)):
          for j in range(1, len(ss)):
              match = self.D[i - 1, j - 1] + (self.match_score if qs[i] == ss[j] else -self.mismatch_penalty)
              self.P[i, j] = max(self.D[i - 1, j] - self.new_gap_penalty - self.continue_gap_penalty, self.P[i - 1, j] - self.continue_gap_penalty)
              self.Q[i, j] = max(self.D[i, j - 1] - self.new_gap_penalty - self.continue_gap_penalty, self.Q[i, j - 1] - self.continue_gap_penalty)
              self.D[i, j] = max(match, self.P[i, j], self.Q[i, j])
              if self.D[i, j] == match: #matrix for traceback based on results from scoring matrix
                  self.pointer[i, j] += 2
              if self.D[i, j] == self.P[i, j]:
                  self.pointer[i, j] += 3
              if self.D[i, j] == self.Q[i, j]:
                  self.pointer[i, j] += 4

        return self.D, self.P, self.Q, self.pointer

    def distance(self, querySequence: str, subjectSequence: str) -> float:
        return super().distance(querySequence, subjectSequence)

    def similarity(self, querySequence: str, subjectSequence: str)->float:
        if querySequence == subjectSequence == "":
            return self.match_score
        D,_, _, _ = self(querySequence, subjectSequence)
        return float(D[D.shape[0]-1,D.shape[1]-1])

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_distance(querySequence, subjectSequence)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_similarity(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64], NDArray[float64], NDArray[float64]]:
        D, P, Q, _ = self(querySequence, subjectSequence)
        return D, P, Q

    def align(self, querySequence: str, subjectSequence: str)->str: 
        _, _, _, pointerMatrix = self(querySequence, subjectSequence)

        qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        i, j = len(qs), len(ss)
        queryAlign, subjectAlign = [], []

        while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
          if pointerMatrix[i,j]in [3, 5, 7, 9]:
              #appends gap and accompanying nucleotide, then moves to the cell above
              subjectAlign.append('-')
              queryAlign.append(qs[i-1])
              i -= 1
          elif pointerMatrix[i,j] in [4, 6, 7, 9]:
              #appends gap and accompanying nucleotide, then moves to the cell to the left
              subjectAlign.append(ss[j-1])
              queryAlign.append('-')
              j -= 1
          elif pointerMatrix[i,j] in [2, 5, 6, 9]:
              #appends match/mismatch then moves to the cell diagonally up and to the left
              queryAlign.append(qs[i-1])
              subjectAlign.append(ss[j-1])
              i -= 1
              j -= 1

        queryAlign = "".join(queryAlign[::-1])
        subjectAlign = "".join(subjectAlign[::-1])

        return f"{queryAlign}\n{subjectAlign}"

class Gotoh_Local(_LOCALBASE):
    def __init__(self, match_score=2, mismatch_penalty=1, new_gap_penalty=3, continue_gap_penalty=2):
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.new_gap_penalty = new_gap_penalty
        self.continue_gap_penalty = continue_gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str) -> tuple[NDArray, NDArray, NDArray]:
        """Compute single alignment matrix"""
        # Initialize matrices
        D = numpy.zeros((len(querySequence) + 1, len(subjectSequence) + 1))
        P = numpy.zeros((len(querySequence) + 1, len(subjectSequence) + 1))
        Q = numpy.zeros((len(querySequence) + 1, len(subjectSequence) + 1))

        # Fill matrices
        for i in range(1, len(querySequence) + 1):
            for j in range(1, len(subjectSequence) + 1):
                score = self.match_score if querySequence[i-1].upper() == subjectSequence[j-1].upper() else -self.mismatch_penalty
                P[i,j] = max(D[i-1,j] - self.new_gap_penalty, P[i-1,j] - self.continue_gap_penalty)
                Q[i,j] = max(D[i,j-1] - self.new_gap_penalty, Q[i,j-1] - self.continue_gap_penalty)
                D[i,j] = max(0, D[i-1,j-1] + score, P[i,j], Q[i,j])

        return D, P, Q

    def _init_alignment_matrices(self, rows: int, cols: int) -> tuple[NDArray, NDArray, NDArray]:
        """Initialize three matrices for Gotoh alignment."""
        return (
            numpy.zeros((rows + 1, cols + 1)),  # D matrix
            numpy.zeros((rows + 1, cols + 1)),  # P matrix
            numpy.zeros((rows + 1, cols + 1))   # Q matrix
        )

    def _fill_cell(self, matrices: tuple[NDArray, NDArray, NDArray], i: int, j: int, char1: str, char2: str) -> None:
        """Fill cell in Gotoh matrices."""
        D, P, Q = matrices
        score = self.match_score if char1.upper() == char2.upper() else -self.mismatch_penalty

        P[i,j] = max(D[i-1,j] - self.new_gap_penalty, 
                     P[i-1,j] - self.continue_gap_penalty)
        Q[i,j] = max(D[i,j-1] - self.new_gap_penalty, 
                     Q[i,j-1] - self.continue_gap_penalty)
        D[i,j] = max(0,
                     D[i-1,j-1] + score,
                     P[i,j],
                     Q[i,j])

    def _get_max_scores(self, matrices_A: tuple[NDArray, NDArray, NDArray],
                       matrices_B: tuple[NDArray, NDArray, NDArray],
                       matrices_AB: tuple[NDArray, NDArray, NDArray]) -> tuple[float, float, float]:
        """Get maximum scores from Gotoh matrices."""
        return matrices_A[0].max(), matrices_B[0].max(), matrices_AB[0].max()

    def distance(self, querySequence: str, subjectSequence: str) -> float:
        return super().distance(querySequence, subjectSequence)

    def similarity(self, querySequence: str, subjectSequence: str)->float:
        if not querySequence and not subjectSequence:
            return 1.0
        matrix, _, _  = self(querySequence, subjectSequence)
        return matrix.max()

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_distance(querySequence, subjectSequence)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        return super().normalized_similarity(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64], NDArray[float64], NDArray[float64]]:
        D, P, Q = self(querySequence, subjectSequence)
        return D, P, Q

    def align(self, querySequence: str, subjectSequence: str)->str:
      matrix, _, _ = self(querySequence, subjectSequence)

      qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
      if matrix.max() == 0:
        return "There is no local alignment!"

      #finds the largest value closest to bottom right of matrix
      i, j = numpy.unravel_index(matrix.argmax(), matrix.shape)

      subjectAlign = []
      queryAlign= []
      score = matrix.max()
      while score > 0:
          score = matrix[i][j]
          if score == 0:
              break
          queryAlign.append(qs[i-1])
          subjectAlign.append(ss[j-1])
          i -= 1
          j -= 1
      queryAlign = "".join(queryAlign[::-1])
      subjectAlign = "".join(subjectAlign[::-1])
      return f"{queryAlign}\n{subjectAlign}"

class Hirschberg():
    def __init__(self, match_score: int = 1, mismatch_penalty: int = 2, gap_penalty: int = 4)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->str:
        qs = ''.join([x.upper() for x in querySequence])
        ss = ''.join([x.upper() for x in subjectSequence])

        if len(qs) == 0:
            return f"{'-' * len(ss)}\n{ss}"
        elif len(ss) == 0:
            return f"{qs}\n{'-' * len(qs)}"
        elif len(qs) == 1 or len(ss) == 1:
            return self._align_simple(qs, ss)

        # Divide and conquer
        xmid = len(qs) // 2

        # Forward score from start to mid
        score_left = self._score(qs[:xmid], ss)
        # Backward score from end to mid
        score_right = self._score(qs[xmid:][::-1], ss[::-1])[::-1]

        # Find optimal split point in subject sequence
        total_scores = score_left + score_right
        ymid = numpy.argmin(total_scores)

        # Recursively align both halves
        left_align = self(qs[:xmid], ss[:ymid])
        right_align = self(qs[xmid:], ss[ymid:])

        # Combine the alignments
        left_q, left_s = left_align.split('\n')
        right_q, right_s = right_align.split('\n')
        return f"{left_q + right_q}\n{left_s + right_s}"

    def _score(self, qs: str, ss: str)->NDArray[float64]:
        # Calculate forward/backward score profile
        prev_row = numpy.zeros(len(ss) + 1, dtype=float64)
        curr_row = numpy.zeros(len(ss) + 1, dtype=float64)

        # Initialize first row
        for j in range(1, len(ss) + 1):
            prev_row[j] = prev_row[j-1] + self.gap_penalty

        # Fill matrix
        for i in range(1, len(qs) + 1):
            curr_row[0] = prev_row[0] + self.gap_penalty
            for j in range(1, len(ss) + 1):
                match_score = (-self.match_score if qs[i-1] == ss[j-1]
                             else self.mismatch_penalty)
                curr_row[j] = min(
                    prev_row[j-1] + match_score,  # match/mismatch
                    prev_row[j] + self.gap_penalty,  # deletion
                    curr_row[j-1] + self.gap_penalty  # insertion
                )
            prev_row, curr_row = curr_row, prev_row

        return prev_row

    def _align_simple(self, qs: str, ss: str)->str:
        score = numpy.zeros((len(qs) + 1, len(ss) + 1), dtype=float64)
        pointer = numpy.zeros((len(qs) + 1, len(ss) + 1), dtype=float64)

        # Initialize first row and column
        for i in range(1, len(qs) + 1):
            score[i,0] = score[i-1,0] + self.gap_penalty
            pointer[i,0] = 1
        for j in range(1, len(ss) + 1):
            score[0,j] = score[0,j-1] + self.gap_penalty
            pointer[0,j] = 2

        # Fill matrices
        for i in range(1, len(qs) + 1):
            for j in range(1, len(ss) + 1):
                match_score = (-self.match_score if qs[i-1] == ss[j-1]
                             else self.mismatch_penalty)
                diag = score[i-1,j-1] + match_score
                up = score[i-1,j] + self.gap_penalty
                left = score[i,j-1] + self.gap_penalty

                score[i,j] = min(diag, up, left)
                if score[i,j] == diag:
                    pointer[i,j] = 3
                elif score[i,j] == up:
                    pointer[i,j] = 1
                else:
                    pointer[i,j] = 2

        # Traceback
        i, j = len(qs), len(ss)
        query_align, subject_align = [], []

        while i > 0 or j > 0:
            if i > 0 and j > 0 and pointer[i,j] == 3:
                query_align.append(qs[i-1])
                subject_align.append(ss[j-1])
                i -= 1
                j -= 1
            elif i > 0 and pointer[i,j] == 1:
                query_align.append(qs[i-1])
                subject_align.append('-')
                i -= 1
            else:
                query_align.append('-')
                subject_align.append(ss[j-1])
                j -= 1

        return f"{''.join(query_align[::-1])}\n{''.join(subject_align[::-1])}"

    def distance(self, querySequence: str, subjectSequence: str)->float:
        """Calculate edit distance between sequences"""
        if not querySequence and not subjectSequence:
            return 0.0
        if not querySequence:
            return self.gap_penalty * len(subjectSequence)
        if not subjectSequence:
            return self.gap_penalty * len(querySequence)

        alignment = self(querySequence, subjectSequence)
        query_align, subject_align = alignment.split('\n')

        dist = 0.0
        for q, s in zip(query_align, subject_align):
            if q == '-' or s == '-':
                dist += self.gap_penalty
            elif q != s:
                dist += self.mismatch_penalty
            # No reduction for matches in distance calculation
        return float(dist)

    def similarity(self, querySequence: str, subjectSequence: str)->float:
        """Calculate similarity score between sequences"""
        if not querySequence and not subjectSequence:
            return 1.0
        if not querySequence or not subjectSequence:
            return 0.0
        alignment = self(querySequence, subjectSequence)
        query_align, subject_align = alignment.split('\n')

        score = 0.0
        for q, s in zip(query_align, subject_align):
            if q == s and q != '-':
                score += self.match_score
            elif q == '-' or s == '-':
                score -= self.gap_penalty
            else:
                score -= self.mismatch_penalty
        return max(0.0, float(score))

    def normalized_distance(self, querySequence: str, subjectSequence: str)->float:
        """Calculate normalized distance between sequences"""
        if not querySequence or not subjectSequence:
            return 1.0
        if querySequence == subjectSequence:
            return 0.0

        raw_dist = self.distance(querySequence, subjectSequence)
        max_len = max(len(querySequence), len(subjectSequence))
        worst_score = max_len * self.mismatch_penalty

        if worst_score == 0:
            return 0.0
        return min(1.0, raw_dist / worst_score)

    def normalized_similarity(self, querySequence: str, subjectSequence: str)->float:
        """Calculate normalized similarity between sequences"""
        if not querySequence or not subjectSequence:
            return 0.0
        if querySequence == subjectSequence:
            return 1.0

        return 1.0 - self.normalized_distance(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->NDArray[float64]:
        if len(querySequence) <= 1 or len(subjectSequence) <= 1:
            score = numpy.zeros((len(querySequence) + 1, len(subjectSequence) + 1), dtype=float64)
            for i in range(len(querySequence) + 1):
                score[i,0] = i * self.gap_penalty
            for j in range(len(subjectSequence) + 1):
                score[0,j] = j * self.gap_penalty
            for i in range(1, len(querySequence) + 1):
                for j in range(1, len(subjectSequence) + 1):
                    match_score = (-self.match_score if querySequence[i-1] == subjectSequence[j-1]
                                 else self.mismatch_penalty)
                    score[i,j] = min(
                        score[i-1,j-1] + match_score,
                        score[i-1,j] + self.gap_penalty,
                        score[i,j-1] + self.gap_penalty
                    )
            return score
        return numpy.array([[]], dtype=float64)

    def align(self, querySequence: str, subjectSequence: str)->str:
        return self(querySequence, subjectSequence)

class Jaro():
    def __init__(self)->None:
        self.match_score = 1
        self.winkler = False
        self.scaling_factor = 1

    def __call__(self, querySequence: str, subjectSequence: str) -> tuple[int, int]:
        qs, ss = (x.upper() for x in [querySequence, subjectSequence])
        if qs == ss:
            return -1, 0
        len1, len2 = len(querySequence), len(subjectSequence)
        max_dist = max(len1, len2)//2 - 1

        matches = 0
        array_qs = [False] * len1
        array_ss = [False] * len2
        for i in range(len1):
            start = max(0, i - max_dist)
            end = min(len2, i + max_dist + 1)
            for j in range(start, end):
                if qs[i] == ss[j] and array_ss[j] == 0:
                    array_qs[i] = array_ss[j] = True
                    matches += 1
                    break
        if matches == 0:
            return 0, 0

        transpositions = 0
        comparison = 0
        for i in range(len1):
            if array_qs[i]:
                while not array_ss[comparison]:
                    comparison += 1
                if qs[i] != ss[comparison]:
                    transpositions += 1
                comparison += 1
        return matches, transpositions//2

    def distance(self, querySequence: str, subjectSequence: str) -> float:
        return 1 - self.similarity(querySequence, subjectSequence)

    def similarity(self, querySequence: str, subjectSequence: str) -> float:
        if not querySequence or not subjectSequence:
            return 1.0 if querySequence == subjectSequence else 0.0

        matches, t = self(querySequence, subjectSequence)
        if matches == 0:
            return 0.0
        if matches == -1:
            return 1.0

        lenqs, lenss = len(querySequence), len(subjectSequence)
        jaro_sim = (1/3)*((matches/lenqs)+(matches/lenss)+((matches-t)/matches))

        if not self.winkler:
            return jaro_sim

        prefix_matches = 0
        max_prefix = min(4, min(lenqs, lenss))
        for i in range(max_prefix):
            if querySequence[i] != subjectSequence[i] or i > len(subjectSequence) - 1:
                break
            prefix_matches += 1
        return jaro_sim + prefix_matches*self.scaling_factor*(1-jaro_sim)

    def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
        return self.distance(querySequence, subjectSequence)

    def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
        return self.similarity(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str) -> NDArray[float64]:
        # dynamic programming variant to show all matches
        qs, ss = [""], [""]
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])
        max_match_dist = max(0, (max(len(ss)-1, len(qs)-1)//2)-1)

        # matrix initialization
        self.alignment_score = numpy.zeros((len(qs), len(ss)))
        for i, query_char in enumerate(qs):
            for j, subject_char in enumerate(ss):
                if i == 0 or j == 0:
                    # keeps first row and column consistent throughout all calculations
                    continue
                dmatch = self.alignment_score[i-1][j-1]
                start = max(1, i - max_match_dist)
                trans_match = ss[start:start + (2 * max_match_dist)]
                if query_char == subject_char or query_char in trans_match:
                    dmatch += 1

                self.alignment_score[i][j] = dmatch
        return self.alignment_score

    def align(self, querySequence: str, subjectSequence: str) -> str:
        """Return aligned sequences showing matches."""
        qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        if qs == ss:
            return f"{''.join(qs)}\n{''.join(ss)}"

        # Initialize arrays for tracking matches
        array_qs = [False] * len(qs)
        array_ss = [False] * len(ss)
        max_dist = max(len(qs), len(ss)) // 2 - 1

        # First pass: mark matches
        for i in range(len(qs)):
            start = max(0, i - max_dist)
            end = min(len(ss), i + max_dist + 1)
            for j in range(start, end):
                if qs[i] == ss[j] and not array_ss[j]:
                    array_qs[i] = array_ss[j] = True
                    break

        # Build global alignment
        queryAlign, subjectAlign = [], []
        i = j = 0

        while i < len(qs) or j < len(ss):
            if i < len(qs) and j < len(ss) and array_qs[i] and array_ss[j] and qs[i] == ss[j]:
                # Add match
                queryAlign.append(qs[i])
                subjectAlign.append(ss[j])
                i += 1
                j += 1
            elif i < len(qs) and not array_qs[i]:
                # Add unmatched query character
                queryAlign.append(qs[i])
                subjectAlign.append('-')
                i += 1
            elif j < len(ss) and not array_ss[j]:
                # Add unmatched subject character
                queryAlign.append('-')
                subjectAlign.append(ss[j])
                j += 1
            elif i < len(qs) and j < len(ss):
                queryAlign.append(qs[i])
                subjectAlign.append(ss[j])
                i += 1
                j += 1
            elif i < len(qs):  # Remaining query characters
                queryAlign.append(qs[i])
                subjectAlign.append('-')
                i += 1
            elif j < len(ss):  # Remaining subject characters
                queryAlign.append('-')
                subjectAlign.append(ss[j])
                j += 1

        return f"{''.join(queryAlign)}\n{''.join(subjectAlign)}"

class Jaro_Winkler(Jaro):
    def __init__(self, scaling_factor = 0.1):
        self.match_score = 1
        self.winkler = True
        #p should not exceed 0.25 else similarity could be larger than 1
        self.scaling_factor = scaling_factor

class Smith_Waterman():
    def __init__(self, match_score:int = 1, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)-> NDArray[float64]:
        qs,ss = [""], [""]
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
            if j == 0 or i == 0:
                #keeps first row and column consistent throughout all calculations
                continue
            if query_char == subject_char:
                match = self.alignment_score[i-1][j-1] + self.match_score
            else:
                match = self.alignment_score[i-1][j-1] - self.mismatch_penalty
            ugap = self.alignment_score[i-1][j] - self.gap_penalty
            lgap = self.alignment_score[i][j-1] - self.gap_penalty
            tmax = max(0, match, lgap, ugap)
            self.alignment_score[i][j] = tmax
        return self.alignment_score

    def distance(self, querySequence: str, subjectSequence: str)->float:
        if not querySequence and not subjectSequence:
            return 0
        return max(map(len, [querySequence, subjectSequence])) - self.similarity(querySequence, subjectSequence)

    def similarity(self, querySequence: str, subjectSequence: str)->float:
        if not querySequence and not subjectSequence:
            return 1.0
        matrix  = self(querySequence, subjectSequence)
        return matrix.max()

    def normalized_distance(self, querySequence: str, subjectSequence: str)->float:
        if not querySequence and not subjectSequence:
            return 0
        dist = self.distance(querySequence, subjectSequence)
        return dist/max(map(len, [querySequence, subjectSequence]))

    def normalized_similarity(self, querySequence: str, subjectSequence: str)->float:
        if not querySequence and not subjectSequence:
            return 1
        similarity = self.similarity(querySequence, subjectSequence)
        return similarity/max(map(len, [querySequence, subjectSequence]))

    def matrix(self, querySequence: str, subjectSequence: str)->NDArray[float64]:
        matrix = self(querySequence, subjectSequence)
        return matrix

    def align(self, querySequence: str, subjectSequence: str)->str:
        matrix = self(querySequence, subjectSequence)

        qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]
        if matrix.max() == 0:
            return "There is no local alignment!"

        #finds the largest value closest to bottom right of matrix
        i, j = list(numpy.where(matrix == matrix.max()))
        i, j = i[-1], j[-1]

        subjectAlign = []
        queryAlign= []
        score = matrix.max()
        while score > 0:
            score = matrix[i][j]
            if score == 0:
                break
            queryAlign.append(qs[i-1])
            subjectAlign.append(ss[j-1])
            i -= 1
            j -= 1
        queryAlign = "".join(queryAlign[::-1])
        subjectAlign = "".join(subjectAlign[::-1])
        return f"{queryAlign}\n{subjectAlign}"

class Longest_Common_Subsequence(_LOCALBASE):
    def __init__(self):
        self.match_score = 1

    def __call__(self, querySequence: str, subjectSequence: str)-> NDArray[float64]:
        qs,ss = [""], [""]
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
            if j == 0 or i == 0:
                #keeps first row and column consistent throughout all calculations
                continue
            if query_char == subject_char:
                match = self.alignment_score[i-1][j-1] + self.match_score
            else:
                match = max(self.alignment_score[i][j-1], self.alignment_score[i-1][j]) 
            self.alignment_score[i][j] = match

        return self.alignment_score

    def distance(self, querysequence: str, subjectsequence: str) -> float:
        return super().distance(querysequence, subjectsequence)

    def similarity(self, querysequence: str, subjectsequence: str) -> float:
        return super().similarity(querysequence, subjectsequence)

    def normalized_distance(self, querysequence: str, subjectsequence: str) -> float:
        return super().normalized_distance(querysequence, subjectsequence)

    def normalized_similarity(self, querysequence: str, subjectsequence: str) -> float:
        return super().normalized_similarity(querysequence, subjectsequence)

    def matrix(self, querySequence: str, subjectSequence: str) -> NDArray:
        return super().matrix(querySequence, subjectSequence)

    def align(self, querySequence: str, subjectSequence: str)->str:
      matrix = self(querySequence, subjectSequence)

      qs = [x.upper() for x in querySequence]
      ss = [x.upper() for x in subjectSequence]
      if matrix.max() == 0:
        return "There is no common subsequence!"

      i, j = len(querySequence), len(subjectSequence)
      common_sub_align = []
      while matrix[i, j] > 0:
          if i == 0 and j == 0:
              break
          if qs[i-1] == ss[j-1]:
              common_sub_align.append(qs[i-1])
              i -= 1
              j -= 1
          elif matrix[i-1,j] >= matrix[i, j-1]:
              i -= 1
          elif matrix[i, j-1] >= matrix[i-1, j]:
              j -= 1
      common_sub_align = "".join(common_sub_align[::-1])
      return f"{common_sub_align}"

class Shortest_Common_Supersequence():
    def __call__(self, querySequence: str, subjectSequence: str)->NDArray[float64]:
        qs,ss = [""], [""]
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        # Matrix initialization with correct shape
        self.alignment_score = numpy.zeros((len(qs),len(ss)), dtype=float64)

        # Fill first row and column
        self.alignment_score[:,0] = [i for i in range(len(qs))]
        self.alignment_score[0,:] = [j for j in range(len(ss))]
        # Fill rest of matrix
        for i in range(1, len(qs)):
            for j in range(1, len(ss)):
                if qs[i] == ss[j]:
                    self.alignment_score[i,j] = self.alignment_score[i-1,j-1]
                else:
                    self.alignment_score[i,j] = min(
                        self.alignment_score[i-1,j] + 1,
                        self.alignment_score[i,j-1] + 1
                    )
        return self.alignment_score

    def distance(self, querySequence: str, subjectSequence: str)->float:
        """Return length of SCS minus length of longer sequence"""
        if not querySequence or not subjectSequence:
            return max(len(querySequence), len(subjectSequence))

        matrix = self(querySequence, subjectSequence)
        return matrix[matrix.shape[0]-1,matrix.shape[1]-1]

    def similarity(self, querySequence: str, subjectSequence: str) -> float:
        """Calculate similarity based on matching positions in supersequence.

        Similarity is the number of positions where characters match between
        the query sequence and the shortest common supersequence.
        """
        if not querySequence or not subjectSequence:
            return 0.0

        scs = self.align(querySequence, subjectSequence)
        return len(scs) - self.distance(querySequence, subjectSequence)

    def normalized_distance(self, querySequence: str, subjectSequence: str)->float:
        """Calculate normalized distance between sequences"""
        if not querySequence or not subjectSequence:
            return 1.0 if (querySequence or subjectSequence) else 0.0
        if querySequence == subjectSequence == "":
            return 0.0
        alignment_len = len(self.align(querySequence, subjectSequence))
        distance = self.distance(querySequence, subjectSequence)
        return distance/alignment_len

    def normalized_similarity(self, querySequence: str, subjectSequence: str)->float:
        """Calculate normalized similarity between sequences"""
        return 1.0 - self.normalized_distance(querySequence, subjectSequence)

    def matrix(self, querySequence: str, subjectSequence: str)->NDArray[float64]:
        return self(querySequence, subjectSequence)

    def align(self, querySequence: str, subjectSequence: str)->str:
        if not querySequence:
            return subjectSequence
        if not subjectSequence:
            return querySequence

        matrix = self(querySequence, subjectSequence)
        qs = [x.upper() for x in querySequence]
        ss = [x.upper() for x in subjectSequence]

        i, j = len(qs), len(ss)
        result = []

        while i > 0 and j > 0:
            if qs[i-1] == ss[j-1]:
                result.append(qs[i-1])
                i -= 1
                j -= 1
            elif matrix[i,j-1] <= matrix[i-1,j]:
                result.append(ss[j-1])
                j -= 1
            else:
                result.append(qs[i-1])
                i -= 1

        # Add remaining characters
        while i > 0:
            result.append(qs[i-1])
            i -= 1
        while j > 0:
            result.append(ss[j-1])
            j -= 1

        return "".join(reversed(result))

hamming = Hamming()
wagner_fischer = Wagner_Fischer()
needleman_wunsch = Needleman_Wunsch()
waterman_smith_beyer = Waterman_Smith_Beyer()
smith_waterman = Smith_Waterman()
hirschberg = Hirschberg()
jaro = Jaro()
jaro_winkler = Jaro_Winkler()
lowrance_wagner = Lowrance_Wagner()
longest_common_subsequence = Longest_Common_Subsequence()
shortest_common_supersequence = Shortest_Common_Supersequence()
gotoh = Gotoh()
gotoh_local = Gotoh_Local()

if __name__ == "__main__":
    main()
