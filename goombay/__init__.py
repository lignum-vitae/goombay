# Base classes
from goombay.algorithms.base import GLOBALBASE, LOCALBASE

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
Jaro_Winkler = editdistance.Jaro_Winkler
Hirschberg = editdistance.Hirschberg
Lowrance_Wagner = editdistance.Lowrance_Wagner
Needleman_Wunsch = editdistance.Needleman_Wunsch
Smith_Waterman = editdistance.Smith_Waterman
Wagner_Fischer = editdistance.Wagner_Fischer
Waterman_Smith_Beyer = editdistance.Waterman_Smith_Beyer
Longest_Common_Subsequence = editdistance.Longest_Common_Subsequence
Shortest_Common_Supersequence = editdistance.Shortest_Common_Supersequence
Gotoh = editdistance.Gotoh
Gotoh_Local = editdistance.Gotoh_Local
