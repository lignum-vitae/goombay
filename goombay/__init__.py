# Base classes
from goombay.algorithms.base import GlobalBase, LocalBase

# Algorithms module
from goombay.algorithms import editdistance
from goombay.algorithms import msa

# Variables from pairwise file
hamming = editdistance.hamming
jaro = editdistance.jaro
jaro_winkler = editdistance.jaro_winkler
hirschberg = editdistance.hirschberg
lowrance_wagner = editdistance.lowrance_wagner
needleman_wunsch = editdistance.needleman_wunsch
gotoh = editdistance.gotoh
gotoh_local = editdistance.gotoh_local
smith_waterman = editdistance.smith_waterman
wagner_fischer = editdistance.wagner_fischer
waterman_smith_beyer = editdistance.waterman_smith_beyer
longest_common_subsequence = editdistance.longest_common_subsequence
longest_common_substring = editdistance.longest_common_substring
shortest_common_supersequence = editdistance.shortest_common_supersequence
lipns = editdistance.lipns
mlipns = editdistance.mlipns

# Variables from multiple sequence alignment file
longest_common_substring_msa = msa.longest_common_substring_msa

# Classes from pairwise file
Hamming = editdistance.Hamming
Jaro = editdistance.Jaro
JaroWinkler = editdistance.JaroWinkler
Hirschberg = editdistance.Hirschberg
LowranceWagner = editdistance.LowranceWagner
NeedlemanWunsch = editdistance.NeedlemanWunsch
Gotoh = editdistance.Gotoh
GotohLocal = editdistance.GotohLocal
SmithWaterman = editdistance.SmithWaterman
WagnerFischer = editdistance.WagnerFischer
WatermanSmithBeyer = editdistance.WatermanSmithBeyer
LongestCommonSubsequence = editdistance.LongestCommonSubsequence
LongestCommonSubstring = editdistance.LongestCommonSubstring
ShortestCommonSupersequence = editdistance.ShortestCommonSupersequence
LIPNS = editdistance.LIPNS
MLIPNS = editdistance.MLIPNS

# Classes from multiple sequence alignment file
LongestCommonSubstringMSA = msa.LongestCommonSubstringMSA
