try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
except ImportError:
    raise ImportError("Numpy is not installed. Please pip install numpy to continue.")

# internal dependencies
from goombay.align.base import LocalBase as _LocalBase
from goombay.align.edit import hamming

__all__ = [
    "LongestCommonSubsequence",
    "longest_common_subsequence",
    "LongestCommonSubstring",
    "longest_common_substring",
    "ShortestCommonSupersequence",
    "shortest_common_supersequence",
    "LIPNS",
    "lipns",
    "MLIPNS",
    "mlipns",
    "LengthRatio",
    "length_ratio",
    "Hamann",
    "hamann",
    "SimpleMatchingCoefficient",
    "simple_matching_coefficient",
    "Prefix",
    "prefix",
    "Postfix",
    "postfix",
    "RatcliffObershelp",
    "ratcliff_obershelp",
]


def main():
    # query = ["WIKIMEDIA", "GESTALT PATTERN MATCHING"]
    # subject = ["WIKIMANIA", "GESTALT PRACTICE"]
    """
    for qs, ss in zip(query, subject):
        print(ratcliff_obershelp.align(qs, ss))
        print(ratcliff_obershelp.similarity(qs, ss))
    """
    query = "HUMAN"
    subject = "CHIMPANZEE"
    print(ratcliff_obershelp(query, subject))


class LongestCommonSubsequence(_LocalBase):
    def __init__(self):
        self.match = 1

    def __call__(self, query_seq: str, subject_seq: str) -> NDArray[float64]:
        qs, ss = [""], [""]
        qs.extend([x.upper() for x in query_seq])
        ss.extend([x.upper() for x in subject_seq])
        qs_len = len(qs)
        ss_len = len(ss)

        # matrix initialisation
        self.score = numpy.zeros((qs_len, ss_len))
        for i in range(1, qs_len):
            for j in range(1, ss_len):
                if qs[i] == ss[j]:
                    match = self.score[i - 1][j - 1] + self.match
                else:
                    match = max(self.score[i][j - 1], self.score[i - 1][j])
                self.score[i][j] = match

        return self.score

    def distance(self, query_seq: str, subject_seq: str) -> float:
        return super().distance(query_seq, subject_seq)

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        return super().similarity(query_seq, subject_seq)

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return super().normalized_distance(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        return super().normalized_similarity(query_seq, subject_seq)

    def matrix(self, query_seq: str, subject_seq: str) -> NDArray:
        return super().matrix(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str) -> list[str]:
        matrix = self(query_seq, subject_seq)

        qs = [x.upper() for x in query_seq]
        ss = [x.upper() for x in subject_seq]

        longest_match = numpy.max(matrix)
        if longest_match <= 1:
            return []

        longest_subseqs = set()
        positions = numpy.argwhere(matrix == longest_match)
        for position in positions:
            temp = []
            i, j = position
            while i != 0 and j != 0:
                if qs[i - 1] == ss[j - 1]:
                    temp.append(qs[i - 1])
                    i -= 1
                    j -= 1
                elif matrix[i - 1, j] >= matrix[i, j - 1]:
                    i -= 1
                elif matrix[i, j - 1] >= matrix[i - 1, j]:
                    j -= 1
            longest_subseqs.add("".join(temp[::-1]))
        return list(longest_subseqs)


class LongestCommonSubstring(_LocalBase):
    def __init__(self):
        self.match = 1

    def __call__(self, query_seq: str, subject_seq: str):
        qs, ss = [""], [""]
        qs.extend([x.upper() for x in query_seq])
        ss.extend([x.upper() for x in subject_seq])
        qs_len = len(qs)
        ss_len = len(ss)

        # matrix initialisation
        alignment_matrix = numpy.zeros((qs_len, ss_len))
        for i in range(1, qs_len):
            for j in range(1, ss_len):
                if qs[i] == ss[j]:
                    match = alignment_matrix[i - 1][j - 1] + self.match
                else:
                    match = 0
                alignment_matrix[i][j] = match
        return alignment_matrix

    def distance(self, query_seq: str, subject_seq: str) -> float:
        return super().distance(query_seq, subject_seq)

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        return super().similarity(query_seq, subject_seq)

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return super().normalized_distance(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        return super().normalized_similarity(query_seq, subject_seq)

    def matrix(self, query_seq: str, subject_seq: str) -> NDArray:
        return super().matrix(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str, min_match: int = 2) -> list[str]:
        matrix = self(query_seq, subject_seq)

        longest_match = numpy.max(matrix)
        if longest_match < min_match or longest_match == 0 or min_match <= 0:
            return [""]

        longest_substrings = []
        positions = numpy.argwhere(matrix == longest_match)
        for position in positions:
            temp = []
            i, j = position
            while matrix[i][j] != 0:
                temp.append(query_seq[i - 1])
                i -= 1
                j -= 1
            longest_substrings.append("".join(temp[::-1]))
        return longest_substrings


class ShortestCommonSupersequence:
    def __call__(self, query_seq: str, subject_seq: str) -> NDArray[float64]:
        qs, ss = [""], [""]
        qs.extend([x.upper() for x in query_seq])
        ss.extend([x.upper() for x in subject_seq])
        qs_len = len(qs)
        ss_len = len(ss)

        # Matrix initialization with correct shape
        self.score = numpy.zeros((qs_len, ss_len), dtype=float64)

        # Fill first row and column
        self.score[:, 0] = [i for i in range(qs_len)]
        self.score[0, :] = [j for j in range(ss_len)]
        # Fill rest of matrix
        for i in range(1, qs_len):
            for j in range(1, ss_len):
                if qs[i] == ss[j]:
                    self.score[i, j] = self.score[i - 1, j - 1]
                else:
                    self.score[i, j] = min(
                        self.score[i - 1, j] + 1,
                        self.score[i, j - 1] + 1,
                    )
        return self.score

    def distance(self, query_seq: str, subject_seq: str) -> float:
        """Return length of SCS minus length of longer sequence"""
        if not query_seq or not subject_seq:
            return max(len(query_seq), len(subject_seq))

        matrix = self(query_seq, subject_seq)
        return matrix[matrix.shape[0] - 1, matrix.shape[1] - 1]

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        """Calculate similarity based on matching positions in supersequence.

        Similarity is the number of positions where characters match between
        the query sequence and the shortest common supersequence.
        """
        if not query_seq or not subject_seq:
            return 0.0

        scs = self.align(query_seq, subject_seq)
        return len(scs) - self.distance(query_seq, subject_seq)

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        """Calculate normalized distance between sequences"""
        if not query_seq or not subject_seq:
            return 1.0 if (query_seq or subject_seq) else 0.0
        if query_seq == subject_seq == "":
            return 0.0
        alignment_len = len(self.align(query_seq, subject_seq))
        distance = self.distance(query_seq, subject_seq)
        return distance / alignment_len

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        """Calculate normalized similarity between sequences"""
        return 1.0 - self.normalized_distance(query_seq, subject_seq)

    def matrix(self, query_seq: str, subject_seq: str) -> NDArray[float64]:
        return self(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str) -> str:
        if not query_seq:
            return subject_seq
        if not subject_seq:
            return query_seq

        matrix = self(query_seq, subject_seq)
        qs = [x.upper() for x in query_seq]
        ss = [x.upper() for x in subject_seq]

        i, j = len(qs), len(ss)
        result = []

        while i > 0 and j > 0:
            if qs[i - 1] == ss[j - 1]:
                result.append(qs[i - 1])
                i -= 1
                j -= 1
            elif matrix[i, j - 1] <= matrix[i - 1, j]:
                result.append(ss[j - 1])
                j -= 1
            else:
                result.append(qs[i - 1])
                i -= 1

        # Add remaining characters
        while i > 0:
            result.append(qs[i - 1])
            i -= 1
        while j > 0:
            result.append(ss[j - 1])
            j -= 1

        return "".join(reversed(result))


class LIPNS:
    # Language-Independent Product Name Search
    def __init__(self, threshold: float = 0.25):
        self.match = 1
        self.threshold = threshold

    def __call__(self, query_seq: str, subject_seq: str):
        qs, ss = [], []
        qs.extend([x.upper() for x in query_seq])
        ss.extend([x.upper() for x in subject_seq])
        qs_len = len(qs)
        ss_len = len(ss)

        # Matrix initialization with correct shape
        score = numpy.zeros((qs_len, ss_len), dtype=float64)

        for i in range(min(qs_len, ss_len)):
            if qs[i] == ss[i]:
                score[i, i] = self.match
        return score

    def distance(self, query_seq: str, subject_seq: str) -> float:
        sim = self.similarity(query_seq, subject_seq)
        return 1 - sim

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq and not subject_seq:
            return 0
        matrix = self(query_seq, subject_seq)
        sim = numpy.sum(matrix)
        sim_score = 1 - (sim / max(len(query_seq), len(subject_seq)))
        return sim_score

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return self.distance(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        return self.similarity(query_seq, subject_seq)

    def matrix(self, query_seq: str, subject_seq: str):
        return self(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str) -> str:
        return f"{query_seq}\n{subject_seq}"

    def is_similar(self, query_seq: str, subject_seq: str) -> int:
        sim_score = self.similarity(query_seq, subject_seq)
        return sim_score <= self.threshold


class MLIPNS(LIPNS):
    # Modified Language-Independent Product Name Search
    def __init__(self, max_mismatch: int = 2):
        self.match = 1
        self.max_mismatch = max_mismatch

    def __call__(self, query_seq: str, subject_seq: str):
        qs, ss = [], []
        qs.extend([x.upper() for x in query_seq])
        ss.extend([x.upper() for x in subject_seq])
        qs_len = len(qs)
        ss_len = len(ss)

        # Matrix initialization with correct shape
        score = numpy.zeros((qs_len, ss_len), dtype=float64)

        if abs(qs_len - ss_len) > self.max_mismatch:
            return score

        i = 0
        max_len = min(qs_len, ss_len)
        mismatch = 0

        while i < max_len:
            if qs[i] != ss[i]:
                mismatch += 1
                if mismatch > self.max_mismatch:
                    return numpy.zeros((qs_len, ss_len), dtype=float64)
            score[i, i] = self.match if qs[i] == ss[i] else 0
            i += 1

        return score

    def distance(self, query_seq: str, subject_seq: str) -> float:
        return super().distance(query_seq, subject_seq)

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        return super().similarity(query_seq, subject_seq)

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return super().normalized_distance(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        return super().normalized_similarity(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str) -> str:
        matrix = self(query_seq, subject_seq)
        if numpy.sum(matrix) == 0:
            return "\n"
        aligned = []
        for i in range(len(query_seq)):
            if matrix[i, i] == 1:
                aligned.append(query_seq[i])
        aligned = "".join(aligned)
        return f"{aligned}\n{aligned}"

    def is_similar(self, query_seq: str, subject_seq: str) -> int:
        if not query_seq and not subject_seq:
            return True
        matrix = self(query_seq, subject_seq)
        if numpy.sum(matrix) > 0:
            return True
        return False


class LengthRatio:
    def __call__(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq and not subject_seq:
            return 1
        if not query_seq or not subject_seq:
            return 0
        query_len = len(query_seq)
        subject_len = len(subject_seq)

        if query_len == subject_len:
            return 1

        ratio = min(query_len, subject_len) / max(query_len, subject_len)
        return ratio

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        return self(query_seq, subject_seq)

    def distance(self, query_seq: str, subject_seq: str) -> float:
        return 1 - self(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        return self(query_seq, subject_seq)

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return 1 - self(query_seq, subject_seq)

    def matrix(self, query_seq: str, subject_seq: str):
        query_len = len(query_seq)
        subject_len = len(subject_seq)
        matrix = numpy.zeros((query_len, subject_len))
        for i in range(query_len):
            for j in range(subject_len):
                if query_seq[i].upper() == subject_seq[j].upper():
                    matrix[i, j] = 1
        return matrix

    def align(self, query_seq: str, subject_seq: str) -> str:
        return f"{query_seq}\n{subject_seq}"


class Hamann:
    def _check_inputs(self, query_seq: str, subject_seq: str) -> None:
        if not isinstance(query_seq, (str)) or not isinstance(subject_seq, (str)):
            raise TypeError("Sequences must be strings")
        if len(query_seq) != len(subject_seq):
            raise ValueError("Sequences must be of equal length")

    def __call__(self, query_seq: str, subject_seq: str, binary: bool):
        self._check_inputs(query_seq, subject_seq)
        if binary:
            matrix = numpy.zeros((2, 2))
            for i in range(len(query_seq)):
                if query_seq[i] == "1" and subject_seq[i] == "1":
                    matrix[0, 0] += 1
                elif query_seq[i] == "1" and subject_seq[i] == "0":
                    matrix[1, 0] += 1
                elif query_seq[i] == "0" and subject_seq[i] == "1":
                    matrix[0, 1] += 1
                elif query_seq[i] == "0" and subject_seq[i] == "0":
                    matrix[1, 1] += 1
        else:
            query_seq = query_seq.upper()
            subject_seq = subject_seq.upper()
            matrix = numpy.zeros((2))
            for i in range(len(query_seq)):
                if query_seq[i] == subject_seq[i]:
                    matrix[0] += 1
                else:
                    matrix[1] += 1
        return matrix

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq or not subject_seq:
            raise ValueError("Sequences must be non-empty")
        matching = hamming.similarity(query_seq, subject_seq)
        non_matching = len(query_seq) - matching
        return (matching - non_matching) / len(query_seq)

    def distance(self, query_seq: str, subject_seq: str) -> float:
        sim = self.similarity(query_seq, subject_seq)
        return 1 - ((sim + 1) / 2)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        sim = self.similarity(query_seq, subject_seq)
        return (sim + 1) / 2

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return self.distance(query_seq, subject_seq)

    def matrix(self, query_seq: str, subject_seq: str, binary: bool):
        return self(query_seq, subject_seq, binary)

    def align(self, query_seq: str, subject_seq: str) -> str:
        return f"{query_seq}\n{subject_seq}"


class SimpleMatchingCoefficient:
    def __call__(self, query_seq: str, subject_seq: str):
        query_seq = query_seq.upper()
        subject_seq = subject_seq.upper()
        matrix = numpy.zeros((1, 1))
        for i in range(len(query_seq)):
            if query_seq[i] == subject_seq[i]:
                matrix[0, 0] += 1
            else:
                matrix[0, 1] += 1
        return matrix

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq or not subject_seq:
            raise ValueError("Strings can not be empty")
        sim = hamming.similarity(query_seq, subject_seq)
        return sim / len(query_seq)

    def distance(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq or not subject_seq:
            raise ValueError("Strings can not be empty")
        dist = hamming.distance(query_seq, subject_seq)
        return dist / len(query_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        return self.similarity(query_seq, subject_seq)

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return self.distance(query_seq, subject_seq)

    def matrix(self, query_seq: str, subject_seq: str):
        return self(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str) -> str:
        return f"{query_seq}\n{subject_seq}"


class Prefix:
    def __call__(self, query_seq: str, subject_seq: str):
        query_seq = query_seq.upper()
        subject_seq = subject_seq.upper()

        query_len = len(query_seq)
        subject_len = len(subject_seq)
        matrix = numpy.zeros((query_len, subject_len))
        for i in range(min(query_len, subject_len)):
            if query_seq[i] != subject_seq[i]:
                break
            matrix[i, i] = 1
        return matrix

    def similarity(self, query_seq: str, subject_seq: str) -> int:
        query_seq = query_seq.upper()
        subject_seq = subject_seq.upper()
        sim = 0
        for i in range(min(len(query_seq), len(subject_seq))):
            if query_seq[i] != subject_seq[i]:
                break
            sim += 1
        return sim

    def distance(self, query_seq: str, subject_seq: str) -> int:
        max_length = max(len(query_seq), len(subject_seq))
        return max_length - self.similarity(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq or not subject_seq:
            raise ValueError("Both strings must be non-empty")
        max_length = max(len(query_seq), len(subject_seq))
        return self.similarity(query_seq, subject_seq) / max_length

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq or not subject_seq:
            raise ValueError("Both strings must be non-empty")
        max_length = max(len(query_seq), len(subject_seq))
        return self.distance(query_seq, subject_seq) / max_length

    def matrix(self, query_seq: str, subject_seq: str):
        return self(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str):
        matrix = self(query_seq, subject_seq)
        alignment = []
        for i in range(min(len(query_seq), len(subject_seq))):
            if matrix[i, i] != 1:
                break
            alignment.append(query_seq[i].upper())
        return "".join(alignment)


class Postfix:
    def __init__(self) -> None:
        self.pre = Prefix()

    def __call__(self, query_seq: str, subject_seq: str):
        query_seq = query_seq.upper()[::-1]
        subject_seq = subject_seq.upper()[::-1]

        query_len = len(query_seq)
        subject_len = len(subject_seq)
        matrix = numpy.zeros((query_len, subject_len))
        for i in range(min(query_len, subject_len)):
            if query_seq[i] != subject_seq[i]:
                break
            matrix[i, i] = 1
        return matrix

    def similarity(self, query_seq: str, subject_seq: str) -> int:
        query_seq = query_seq[::-1]
        subject_seq = subject_seq[::-1]
        return self.pre.similarity(query_seq, subject_seq)

    def distance(self, query_seq: str, subject_seq: str) -> int:
        max_length = max(len(query_seq), len(subject_seq))
        return max_length - self.similarity(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq or not subject_seq:
            raise ValueError("Both strings must be non-empty")
        max_length = max(len(query_seq), len(subject_seq))
        return self.similarity(query_seq, subject_seq) / max_length

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq or not subject_seq:
            raise ValueError("Both strings must be non-empty")
        max_length = max(len(query_seq), len(subject_seq))
        return self.distance(query_seq, subject_seq) / max_length

    def matrix(self, query_seq: str, subject_seq: str):
        return self(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str):
        matrix = self(query_seq, subject_seq)
        alignment = []
        for i in range(min(len(query_seq), len(subject_seq))):
            if matrix[i, i] != 1:
                break
            alignment.append(query_seq[-i - 1].upper())
        return "".join(alignment[::-1])


class RatcliffObershelp:
    def __call__(self, query_seq: str, subject_seq: str):
        matched = []
        stack = [(query_seq.upper(), subject_seq.upper())]
        while len(stack) >= 1:
            # Get remaining characters
            remaining = stack.pop()
            qs = remaining[0]
            ss = remaining[1]

            # Find LSString
            matches = longest_common_substring.align(qs, ss, min_match=1)
            if matches == [""]:
                continue
            # Save matches and add remaining back to stack
            q_idx = qs.find(matches[0])
            s_idx = ss.find(matches[0])
            left = (qs[:q_idx], ss[:s_idx])
            right = (qs[q_idx + len(matches[0]) :], ss[s_idx + len(matches[0]) :])
            stack.extend([left, right])
            matched.append(matches[0])
        return matched

    def distance(self, query_seq: str, subject_seq: str) -> float:
        return 1 - self.similarity(query_seq, subject_seq)

    def similarity(self, query_seq: str, subject_seq: str) -> float:
        if not query_seq and not subject_seq:
            return 1.0
        if not query_seq or not subject_seq:
            return 0.0

        sim = self.align(query_seq, subject_seq)
        matches = sum(len(s) for s in sim)
        return (2 * matches) / (len(query_seq) + len(subject_seq))

    def normalized_distance(self, query_seq: str, subject_seq: str) -> float:
        return self.distance(query_seq, subject_seq)

    def normalized_similarity(self, query_seq: str, subject_seq: str) -> float:
        return self.similarity(query_seq, subject_seq)

    def align(self, query_seq: str, subject_seq: str) -> list[str]:
        forward = self(query_seq, subject_seq)
        reverse = self(subject_seq, query_seq)
        forward_len = sum(len(f) for f in forward)
        reverse_len = sum(len(r) for r in reverse)
        return forward if forward_len > reverse_len else reverse


longest_common_subsequence = LongestCommonSubsequence()
longest_common_substring = LongestCommonSubstring()
shortest_common_supersequence = ShortestCommonSupersequence()
lipns = LIPNS()
mlipns = MLIPNS()
length_ratio = LengthRatio()
hamann = Hamann()
simple_matching_coefficient = SimpleMatchingCoefficient()
prefix = Prefix()
postfix = Postfix()
ratcliff_obershelp = RatcliffObershelp()

if __name__ == "__main__":
    main()
