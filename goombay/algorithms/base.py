class GLOBALBASE():
    def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
      matrix, _ = self(querySequence, subjectSequence)
      return matrix

    def distance(self, querySequence: str, subjectSequence: str)->float:
      matrix, _ = self(querySequence, subjectSequence)
      return matrix[matrix.shape[0]-1,matrix.shape[1]-1]

    def similarity(self, querySequence: str, subjectSequence: str)->float:
      return max(len(querySequence),len(subjectSequence)) - self.distance(querySequence, subjectSequence)

    def normalized_distance(self, querySequence: str, subjectSequence: str)->float:
      dist = self.distance(querySequence, subjectSequence)
      return dist/max(map(len, [querySequence, subjectSequence]))

    def normalized_similarity(self, querySequence: str, subjectSequence: str)->float:
      sim = self.similarity(querySequence, subjectSequence)
      return sim/max(map(len, [querySequence, subjectSequence]))

    def align(self, querySequence: str, subjectSequence: str)->str: 
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

class LOCALBASE():
    def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
      matrix = self(querySequence, subjectSequence)
      return matrix

    def similarity(self, querySequence: str, subjectSequence: str)->float:
      matrix  = self(querySequence, subjectSequence)
      return matrix.max()

    def distance(self, querySequence: str, subjectSequence: str)->float:
      return max(map(len, [querySequence, subjectSequence])) - self.similarity(querySequence, subjectSequence)

    def normalized_distance(self, querySequence: str, subjectSequence: str)->float:
      dist = self.distance(querySequence, subjectSequence)
      return dist/max(map(len, [querySequence, subjectSequence]))

    def normalized_similarity(self, querySequence: str, subjectSequence: str)->float:
      similarity = self.similarity(querySequence, subjectSequence)
      return similarity/max(map(len, [querySequence, subjectSequence]))
