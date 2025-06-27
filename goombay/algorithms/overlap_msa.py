def main():
    seq1 = "HOUSEOFCARDSFALLDOWN"
    seq2 = "HOUSECARDFALLDOWN"
    seq3 = "FALLDOWN"
    print(longest_common_substring_msa.align([seq2, seq1, seq3]))
    print(longest_common_substring_msa.distance([seq2, seq1, seq3]))
    print(longest_common_substring_msa.normalized_distance([seq2, seq1, seq3]))
    print(longest_common_substring_msa.similarity([seq2, seq1, seq3]))
    print(longest_common_substring_msa.normalized_similarity([seq2, seq1, seq3]))


__all__ = ["LongestCommonSubstringMSA", "longest_common_substring_msa"]


class LongestCommonSubstringMSA:
    def _common_substrings(self, seqs: list[str]) -> list[str]:
        if not isinstance(seqs, list):
            raise TypeError("common_substrings expects a list of strings")
        if len(seqs) != 2:
            raise ValueError("common_substrings requires exactly two strings")

        s1, s2 = seqs
        common = set()
        s2_length = len(s2)

        for start in range(s2_length):
            for end in range(start + 2, s2_length + 1):  # substrings of length >= 2
                substr = s2[start:end]
                if substr not in s1:
                    break
                common.add(substr)

        return list(common)

    def __call__(self, seqs: list[str]) -> list[str]:
        if (not isinstance(seqs, list) and not isinstance(seqs, tuple)) or not all(
            isinstance(s, str) for s in seqs
        ):
            raise TypeError("longest_common_substring_msa expects a list of strings")
        if len(seqs) < 2:
            raise ValueError("Provide at least two seqs")

        seqs = [seq.upper() for seq in seqs]

        # Generate substrings from first and last strings
        motifs = self._common_substrings([seqs[0], seqs[-1]])
        if not motifs:
            return [""]

        longest = []
        longest_len = -1
        motifs.sort(key=len, reverse=True)
        for motif in motifs:
            motif_len = len(motif)
            if all(motif in seq for seq in seqs) and motif_len >= longest_len:
                longest.append(motif)
                longest_len = motif_len
            if motif_len < longest_len:
                break
        return longest

    def align(self, seqs: list[str]) -> list[str]:
        return self(seqs)

    def distance(self, seqs: list[str]) -> int:
        max_match = len(max(seqs, key=len))
        if any(not seq for seq in seqs):
            return max_match
        if max_match <= 1 and all(seqs[0] in seq for seq in seqs):
            return 0
        return max_match - self.similarity(seqs)

    def similarity(self, seqs: list[str]) -> int:
        lcsub = self(seqs)
        if not lcsub:
            return 0
        if all(len(seq) == 0 for seq in seqs) or all(seqs[0] == seq for seq in seqs):
            return len(seqs[0]) if len(seqs[0]) >= 1 else 1
        return len(lcsub[0])

    def normalized_distance(self, seqs: list[str]) -> float:
        return 1 - self.normalized_similarity(seqs)

    def normalized_similarity(self, seqs: list[str]) -> float:
        if len(max(seqs, key=len)) == 1 and all(seqs[0] in seq for seq in seqs):
            return 1.0
        lcsub = self(seqs)
        if not lcsub:
            return 0.0
        return len(lcsub[0]) / len(min(seqs, key=len))


longest_common_substring_msa = LongestCommonSubstringMSA()

if __name__ == "__main__":
    main()
