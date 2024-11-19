#built-in
from __future__ import annotations

#internal dependencies
from goombay.algorithms.base import GLOBALBASE as __GLOBALBASE, LOCALBASE as __LOCALBASE

try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
except ImportError:
    raise ImportError("Please pip install all dependencies from requirements.txt!")

def main():
    qqs = "WHORUNSTHEWORLDGIRLS"
    sss = "DOGSSTHEWORLDEXUBERENT"

    print(gotoh_local(qqs, sss)) 
    print(gotoh_local.align(qqs, sss)) 

class Wagner_Fischer(__GLOBALBASE): #Levenshtein Distance
    def __init__(self)->None:
        self.gap_penalty = 1

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

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
              substitution_cost = 0
              if i == 0 or j == 0:
                  #keeps first row and column consistent throughout all calculations
                  continue
              if query_char != subject_char:
                  substitution_cost = 1
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

class Lowrance_Wagner(__GLOBALBASE): #Damerau-Levenshtein distance
    def __init__(self)->None:
        self.gap_penalty = 1

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

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
              substitution_cost = 0
              if i == 0 or j == 0:
                  #keeps first row and column consistent throughout all calculations
                  continue
              if query_char != subject_char:
                  substitution_cost = 1
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
          elif pointerMatrix[i,j] in [8, 10, 11, 12, 13, 14, 15, 17]:
              queryAlign.extend([qs[i-1],qs[i-2]])
              subjectAlign.extend([ss[j-2],ss[j-1]])
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

class Hamming(__GLOBALBASE):
    def __int_pair(self, querySequence: str|int, subjectSequence: str|int) -> bool:
      querySequence, subjectSequence = str(querySequence), str(subjectSequence)
      if querySequence.isalpha() and subjectSequence.isalpha():
          return False
      if querySequence.isdigit() and subjectSequence.isdigit():
          return True 
      raise ValueError("Both sequences must be either all letters or all numbers")
        
    def align(self, querySequence: str|int, subjectSequence: str|int)->str:
        if self.__int_pair(querySequence, subjectSequence):
            qs, ss = int(querySequence), int(subjectSequence)
            return f"{bin(qs)}\n{bin(ss)}"
        return f"{querySequence}\n{subjectSequence}"

    def matrix(self, qs: str, ss: str) -> None:
        return None

    def __call__(self, querySequence: str|int, subjectSequence: str|int)->tuple[int,list[int]]:
        if self.__int_pair(querySequence, subjectSequence):
            qs, ss = bin(querySequence), bin(subjectSequence)
        else:
            qs,ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]

        if len(qs) == 1 and len(ss) == 1:
            dist = 1 if qs != ss else 0
            dist_array = [dist]
            return dist, dist_array

        shortlen = min(map(len, [ss,qs]))
        short = qs if len(qs) == shortlen else ss
        long = ss if len(qs) == shortlen else qs

        dist = 0
        dist_array = []
        for i, char in enumerate(short):
            if char != long[i]:
                dist += 1
                dist_array.append(0)
                continue
            dist_array.append(1)

        dist += len(long)-len(short)
        dist_array.extend([1]*(len(long)-len(short)))
        return dist, dist_array

    def distance(self, querySequence: str|int, subjectSequence: str|int)->int:
        if self.__int_pair(querySequence, subjectSequence):
            qs, ss = int(querySequence), int(subjectSequence)
            return bin(qs ^ ss).count("1")
        query = set([(x, y) for (x, y) in enumerate(querySequence)]) 
        subject = set([(x, y) for (x, y) in enumerate(subjectSequence)]) 
        qs,sq = query-subject, subject-query
        dist = max(map(len,[qs,sq]))
        return dist 

    def binary_distance_array(self, querySequence: str, subjectSequence: str)->list[int]:
        _, distarray = self(querySequence, subjectSequence)
        return distarray

    def binary_similarity_array(self, querySequence: str, subjectSequence: str)->list[int]:
        _, distarray = self(querySequence, subjectSequence)
        simarray = [1 if num == 0 else 0 for num in distarray]
        return simarray

class Needleman_Wunsch(__GLOBALBASE):
    def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
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
        self.alignment_score[:,0] = [n*self.gap_penalty for n in range(len(qs))]
        self.alignment_score[0,:] = [n*self.gap_penalty for n in range(len(ss))]

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
              if i == 0 or j == 0:
                  #keeps first row and column consistent throughout all calculations
                  continue
              if query_char == subject_char:
                  match = self.alignment_score[i-1][j-1] - self.match_score
              else:
                  match = self.alignment_score[i-1][j-1] + self.mismatch_penalty
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

class Waterman_Smith_Beyer(__GLOBALBASE):
    def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, new_gap_penalty:int = 1, continue_gap_penalty:int = 1)->None:
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
        self.alignment_score[:,0] = [self.new_gap_penalty + n * self.continue_gap_penalty for n in range(len(qs))]
        self.alignment_score[0,:] = [self.new_gap_penalty + n * self.continue_gap_penalty for n in range(len(ss))] 
        self.alignment_score[0][0] = 0

        for i, subject in enumerate(qs):
            for j, query in enumerate(ss):
                if i == 0 or j == 0:
                    #keeps first row and column consistent throughout all calculations
                    continue
                if subject == query: 
                    matchScore = self.alignment_score[i-1][j-1] - self.match_score
                else:
                    matchScore = self.alignment_score[i-1][j-1] + self.mismatch_penalty
                #both gaps defaulted to continue gap penalty
                ugapScore = self.alignment_score[i-1][j] + self.continue_gap_penalty
                lgapScore = self.alignment_score[i][j-1] + self.continue_gap_penalty
                #if cell before i-1 or j-1 is gap, then this is a gap continuation
                if self.alignment_score[i-1][j] != (self.alignment_score[i-2][j]) + self.new_gap_penalty + self.continue_gap_penalty:
                    ugapScore += self.new_gap_penalty
                if self.alignment_score[i][j-1] != (self.alignment_score[i][j-2]) + self.new_gap_penalty + self.continue_gap_penalty:
                    lgapScore += self.new_gap_penalty
                tmin = min(matchScore, lgapScore, ugapScore)

                self.alignment_score[i][j] = tmin #lowest value is best choice
                #matrix for traceback based on results from scoring matrix
                if matchScore == tmin: 
                    self.pointer[i,j] += 2
                elif ugapScore == tmin:
                    self.pointer[i,j] += 3
                elif lgapScore == tmin:
                    self.pointer[i,j] += 4
        return self.alignment_score, self.pointer

class Gotoh(__GLOBALBASE):
    def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, new_gap_penalty:int = 2, continue_gap_penalty: int = 1)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.new_gap_penalty = new_gap_penalty
        self.continue_gap_penalty = continue_gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64],NDArray[float64],NDArray[float64]]:
        qs,ss = [""], [""] 
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.D = numpy.full((len(qs),len(ss)), numpy.inf)
        self.P = numpy.full((len(qs), len(ss)), numpy.inf)
        self.P[:,0] = 0
        self.Q = numpy.full((len(qs), len(ss)), numpy.inf)
        self.Q[0,:] = 0
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.D[:,0] = [(self.continue_gap_penalty * n + self.new_gap_penalty) for n in range(len(qs))]
        self.D[0,:] = [(self.continue_gap_penalty * n + self.new_gap_penalty) for n in range(len(ss))]
        self.D[0, 0] = 0

        for i in range(1, len(qs)):
          for j in range(1, len(ss)):
              match = self.D[i - 1, j - 1] + (-self.match_score if qs[i] == ss[j] else self.mismatch_penalty)
              self.P[i, j] = min(self.D[i - 1, j] + self.new_gap_penalty + self.continue_gap_penalty, self.P[i - 1, j] + self.continue_gap_penalty)
              self.Q[i, j] = min(self.D[i, j - 1] + self.new_gap_penalty + self.continue_gap_penalty, self.Q[i, j - 1] + self.continue_gap_penalty)
              self.D[i, j] = min(match, self.P[i, j], self.Q[i, j])
              if self.D[i, j] == match: #matrix for traceback based on results from scoring matrix
                  self.pointer[i, j] += 2
              if self.D[i, j] == self.P[i, j]:
                  self.pointer[i, j] += 3
              if self.D[i, j] == self.Q[i, j]:
                  self.pointer[i, j] += 4
        return self.D, self.P, self.Q, self.pointer

    def matrix(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64], NDArray[float64], NDArray[float64]]:
        D, P, Q, _ = self(querySequence, subjectSequence)
        return D, P, Q

    def distance(self, querySequence: str, subjectSequence: str)->float:
      D,_, _, _ = self(querySequence, subjectSequence)
      return float(D[D.shape[0]-1,D.shape[1]-1])

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

class Gotoh_Local(__LOCALBASE):
    def __init__(self, match_score:int = 1, mismatch_penalty:int = 1, new_gap_penalty:int = 2, continue_gap_penalty: int = 1)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.new_gap_penalty = new_gap_penalty
        self.continue_gap_penalty = continue_gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64],NDArray[float64]]:
        qs,ss = [""], [""] 
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.D = numpy.full((len(qs),len(ss)), -numpy.inf)
        self.P = numpy.full((len(qs), len(ss)), -numpy.inf)
        self.P[:,0] = 0
        self.Q = numpy.full((len(qs), len(ss)), -numpy.inf)
        self.Q[0,:] = 0
        #initialisation of starter values for first column and first row
        self.D[:,0] = 0
        self.D[0,:] = 0

        for i in range(1, len(qs)):
          for j in range(1, len(ss)):
              match = self.D[i - 1, j - 1] + (self.match_score if qs[i] == ss[j] else -self.mismatch_penalty)
              self.P[i, j] = max(self.D[i - 1, j] - self.new_gap_penalty - self.continue_gap_penalty, self.P[i - 1, j] - self.continue_gap_penalty)
              self.Q[i, j] = max(self.D[i, j - 1] - self.new_gap_penalty - self.continue_gap_penalty, self.Q[i, j - 1] - self.continue_gap_penalty)
              self.D[i, j] = max(match, self.P[i, j], self.Q[i, j], 0)
        return self.D, self.P, self.Q

    def matrix(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64], NDArray[float64], NDArray[float64]]:
        D, P, Q = self(querySequence, subjectSequence)
        return D, P, Q

    def align(self, querySequence: str, subjectSequence: str)->str:
      matrix, _, _ = self(querySequence, subjectSequence)

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

    def similarity(self, querySequence: str, subjectSequence: str)->float:
      matrix, _, _  = self(querySequence, subjectSequence)
      return matrix.max()

class Hirschberg():
    def __init__(self, match_score: int = 1, mismatch_penalty: int = -1, gap_penalty: int = -2)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def __call__(self, querySequence: str, subjectSequence:str)->tuple[str,str]:
        qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])

        if len(qs) == 0:
            return '-' * len(ss), ss
        elif len(ss) == 0:
            return qs, '-' * len(qs)
        elif len(qs) == 1 or len(ss) == 1:
            return self._align(qs, ss)
        else:
            xlen = len(qs)
            xmid = xlen // 2
            score_l = self._score(qs[:xmid], ss)
            score_r = self._score(qs[xmid:][::-1], ss[::-1])[::-1]
            ymid = numpy.argmax(score_l + score_r)

            A_left, B_left = self(qs[:xmid], ss[:ymid])
            A_right, B_right = self(qs[xmid:], ss[ymid:])
            return A_left + A_right, B_left + B_right

    def _score(self, qs, ss):
        m, n = len(qs), len(ss)
        score = numpy.zeros((2, n + 1))

        for j in range(1, n + 1):
            score[0][j] = score[0][j - 1] + self.gap_penalty
        for i in range(1, m + 1):
            score[1][0] = score[0][0] + self.gap_penalty
            for j in range(1, n + 1):
                match = score[0][j - 1] + (self.match_score if qs[i - 1] == ss[j - 1] else self.mismatch_penalty)
                delete = score[0][j] + self.gap_penalty
                insert = score[1][j - 1] + self.gap_penalty
                score[1][j] = max(match, delete, insert)
            score[0] = score[1]
        return score[1]

    def _align(self, qs, ss):
        m, n = len(qs), len(ss)
        score = numpy.zeros((m + 1, n + 1))
        pointer = numpy.zeros((m + 1, n + 1))

        for i in range(1, m + 1):
            score[i][0] = score[i - 1][0] + self.gap_penalty
            pointer[i][0] = 1
        for j in range(1, n + 1):
            score[0][j] = score[0][j - 1] + self.gap_penalty
            pointer[0][j] = 2
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if qs[i-1] == ss[j-1]: 
                    match = score[i-1][j-1] + self.match_score
                else:
                    match = score[i-1][j-1] + self.mismatch_penalty
                delete = score[i - 1][j] + self.gap_penalty
                insert = score[i][j - 1] + self.gap_penalty
                score[i][j] = max(match, delete, insert)
                if score[i][j] == match:
                    pointer[i][j] = 3
                elif score[i][j] == delete:
                    pointer[i][j] = 1
                else:
                    pointer[i][j] = 2

        queryAlign, subjectAlign = "", ""
        i, j = m, n
        while i > 0 or j > 0:
            if pointer[i][j] == 3:
                queryAlign = qs[i - 1] + queryAlign
                subjectAlign = ss[j - 1] + subjectAlign
                i -= 1
                j -= 1
            elif pointer[i][j] == 1:
                queryAlign = qs[i - 1] + queryAlign
                subjectAlign = '-' + subjectAlign
                i -= 1
            else:
                queryAlign = '-' + queryAlign
                subjectAlign = ss[j - 1] + subjectAlign
                j -= 1
        return queryAlign, subjectAlign
    def align(self, qs: str, ss: str) -> str:
        queryAlign, subjectAlign = self(qs, ss)
        return f"{queryAlign}\n{subjectAlign}"

class Jaro():
  def __init__(self):
    self.match_score = 1
    self.winkler = False
      
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

  def similarity(self, querySequence: str, subjectSequence: str) -> float:
      matches, t = self(querySequence, subjectSequence)
      if matches == 0:
          return 0.0
      if matches == -1:
          return 1.0
      jaro_sim = (1/3)*((matches/len(querySequence))+(matches/len(subjectSequence))+((matches-t)/matches))
      if not self.winkler:
          return jaro_sim
      prefix_matches = 0
      for i in range(4):
          if querySequence[i] != subjectSequence[i] or i > len(subjectSequence) - 1:
              break
          prefix_matches += 1
      return jaro_sim + prefix_matches*self.scaling_factor*(1-jaro_sim)

  def normalized_similarity(self, querySequence: str, subjectSequence: str) -> float:
      return round(self.similarity(querySequence, subjectSequence), 2)
      
  def distance(self, querySequence: str, subjectSequence: str) -> float:
      return 1 - self.similarity(querySequence, subjectSequence)

  def normalized_distance(self, querySequence: str, subjectSequence: str) -> float:
      return round(self.distance(querySequence, subjectSequence), 2)

  def matrix(self, querySequence: str, subjectSequence: str) -> NDArray[float64]:
    #dynamic programming variant to show all matches
    qs,ss = [""], [""] 
    qs.extend([x.upper() for x in querySequence])
    ss.extend([x.upper() for x in subjectSequence])
    max_match_dist = max(0, (max(len(ss)-1, len(qs)-1)//2)-1)

    #matrix initialisation
    self.alignment_score = numpy.zeros((len(qs),len(ss)))
    for i, query_char in enumerate(qs):
      for j, subject_char in enumerate(ss):
          if i == 0 or j == 0:
              #keeps first row and column consistent throughout all calculations
              continue
          dmatch = self.alignment_score[i-1][j-1]
          start = max(1, i-max_match_dist)
          trans_match = ss[start:start+(2*max_match_dist)]
          if query_char == subject_char or query_char in trans_match:
            dmatch += 1

          self.alignment_score[i][j] = dmatch
    return self.alignment_score

class Jaro_Winkler(Jaro):
    def __init__(self, scaling_factor = 0.1):
        self.match_score = 1
        self.winkler = True
        #p should not exceed 0.25 else similarity could be larger than 1
        self.scaling_factor = scaling_factor

class Smith_Waterman(__LOCALBASE):
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

class Longest_Common_Subsequence(__LOCALBASE):
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
 
    def align(self, querySequence: str, subjectSequence: str)->str:
      matrix = self(querySequence, subjectSequence)

      qs = [x.upper() for x in querySequence]
      if matrix.max() == 0:
        return "There is no common subsequence!"

      i, j = len(querySequence), len(subjectSequence)
      common_sub_align = []
      while matrix[i, j] > 0:
          if i == 0 and j == 0:
              break
          if querySequence[i-1] == subjectSequence[j-1]:
              common_sub_align.append(qs[i-1])
              i -= 1
              j -= 1
          elif matrix[i-1,j] >= matrix[i, j-1]:
              i -= 1
          elif matrix[i, j-1] >= matrix[i-1, j]:
              j -= 1
      common_sub_align = "".join(common_sub_align[::-1])
      return f"{common_sub_align}"

class Shortest_Common_Supersequence(__LOCALBASE):
    def __init__(self):
        self.match_score = 1

    def __call__(self, querySequence: str, subjectSequence: str)-> NDArray[float64]: 
        qs,ss = [""], [""] 
        qs.extend([x.upper() for x in querySequence])
        ss.extend([x.upper() for x in subjectSequence])

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss))) 
        self.alignment_score[:,0] = [x for x in range(len(qs))]
        self.alignment_score[0,:] = [x for x in range(len(ss))]

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
            if j == 0 or i == 0:
                #keeps first row and column consistent throughout all calculations
                continue
            if query_char == subject_char: 
                match = self.alignment_score[i-1][j-1] + self.match_score
            else:
                match = min(self.alignment_score[i][j-1], self.alignment_score[i-1][j]) + self.match_score
            self.alignment_score[i][j] = match
        return self.alignment_score

    def align(self, querySequence: str, subjectSequence: str)->str:
      matrix = self(querySequence, subjectSequence)

      qs, ss = [x.upper() for x in querySequence], [x.upper() for x in subjectSequence]

      i, j = len(querySequence), len(subjectSequence)
      common_super_align= []
      while i * j > 0:
          if querySequence[i-1] == subjectSequence[j-1]:
              common_super_align.append(qs[i-1])
              i -= 1
              j -= 1
          elif matrix[i-1,j] > matrix[i,j-1]:
              common_super_align.append(ss[j-1])
              j -= 1
          else:
              common_super_align.append(qs[i-1])
              i -= 1
      while i > 0:
          common_super_align.append(qs[i-1])
          i -= 1
      while j > 0:
          common_super_align.append(ss[j-1])
          j -= 1

      common_super_align = "".join(common_super_align[::-1])
      return f"{common_super_align}"

    def similarity(self, querySequence: str, subjectSequence: str)->float:
      matrix  = longest_common_subsequence(querySequence, subjectSequence)
      return matrix.max()

    def distance(self, querySequence: str, subjectSequence: str)->float:
      matrix = self(querySequence, subjectSequence)
      sim = self.similarity(querySequence, subjectSequence)
      return (matrix.max() - sim)/2

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
