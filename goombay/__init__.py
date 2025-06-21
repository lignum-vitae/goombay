# Base classes
from goombay.algorithms.base import GlobalBase, LocalBase

# Algorithms module
from goombay.algorithms import edit
from goombay.algorithms import overlap
from goombay.algorithms import overlap_msa

# Variables from edit-based file
hamming = edit.hamming
jaro = edit.jaro
jaro_winkler = edit.jaro_winkler
hirschberg = edit.hirschberg
lowrance_wagner = edit.lowrance_wagner
needleman_wunsch = edit.needleman_wunsch
gotoh = edit.gotoh
gotoh_local = edit.gotoh_local
smith_waterman = edit.smith_waterman
wagner_fischer = edit.wagner_fischer
waterman_smith_beyer = edit.waterman_smith_beyer

# variables from overlap-based file
longest_common_subsequence = overlap.longest_common_subsequence
longest_common_substring = overlap.longest_common_substring
shortest_common_supersequence = overlap.shortest_common_supersequence
lipns = overlap.lipns
mlipns = overlap.mlipns

# Variables from multiple sequence alignment file
longest_common_substring_msa = overlap_msa.longest_common_substring_msa

# Classes from edit-based file
Hamming = edit.Hamming
Jaro = edit.Jaro
JaroWinkler = edit.JaroWinkler
Hirschberg = edit.Hirschberg
LowranceWagner = edit.LowranceWagner
NeedlemanWunsch = edit.NeedlemanWunsch
Gotoh = edit.Gotoh
GotohLocal = edit.GotohLocal
SmithWaterman = edit.SmithWaterman
WagnerFischer = edit.WagnerFischer
WatermanSmithBeyer = edit.WatermanSmithBeyer

# Classes from overlap-based file
LongestCommonSubsequence = overlap.LongestCommonSubsequence
LongestCommonSubstring = overlap.LongestCommonSubstring
ShortestCommonSupersequence = overlap.ShortestCommonSupersequence
LIPNS = overlap.LIPNS
MLIPNS = overlap.MLIPNS

# Classes from multiple sequence alignment file
LongestCommonSubstringMSA = overlap_msa.LongestCommonSubstringMSA
