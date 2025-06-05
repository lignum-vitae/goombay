# Base classes
from goombay.algorithms.base import GlobalBase, LocalBase

# Algorithms package
from goombay.algorithms import editdistance

# Variables from algorithms package
hamming = editdistance.hamming
jaro = editdistance.jaro
jaro_winkler = editdistance.jaro_winkler
hirschberg = editdistance.hirschberg
lowrance_wagner = editdistance.lowrance_wagner
needleman_wunsch = editdistance.needleman_wunsch
smith_waterman = editdistance.smith_waterman
wagner_fischer = editdistance.wagner_fischer
waterman_smith_beyer = editdistance.waterman_smith_beyer
longest_common_subsequence = editdistance.longest_common_subsequence
shortest_common_supersequence = editdistance.shortest_common_supersequence
gotoh = editdistance.gotoh
gotoh_local = editdistance.gotoh_local

# Classes from algorithms package
Hamming = editdistance.Hamming
Jaro = editdistance.Jaro
JaroWinkler = editdistance.JaroWinkler
Hirschberg = editdistance.Hirschberg
LowranceWagner = editdistance.LowranceWagner
NeedlemanWunsch = editdistance.NeedlemanWunsch
SmithWaterman = editdistance.SmithWaterman
WagnerFischer = editdistance.WagnerFischer
WatermanSmithBeyer = editdistance.WatermanSmithBeyer
LongestCommonSubsequence = editdistance.LongestCommonSubsequence
ShortestCommonSupersequence = editdistance.ShortestCommonSupersequence
Gotoh = editdistance.Gotoh
GotohLocal = editdistance.GotohLocal
