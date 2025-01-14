[![Static Badge](https://img.shields.io/badge/Project_Name-Goombay-blue)](https://github.com/dawnandrew100/goombay)
[![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Flignum-vitae%2Fgoombay%2Fmaster%2Fpyproject.toml)](https://github.com/lignum-vitae/goombay/blob/master/pyproject.toml)
[![PyPI version](https://img.shields.io/pypi/v/goombay.svg)](https://pypi.python.org/pypi/goombay)
[![License](https://img.shields.io/pypi/l/goombay.svg)](https://github.com/dawnandrew100/goombay/blob/master/LICENSE)
[![GitHub branch check runs](https://img.shields.io/github/check-runs/dawnandrew100/goombay/master)](https://github.com/dawnandrew100/goombay)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/goombay)](https://pypi.python.org/pypi/goombay)

# Goombay
This Python project contains several sequence alignment algorithms that can also produce scoring matrices for Needleman-Wunsch, Gotoh, Smith-Waterman, Wagner-Fischer, Waterman-Smith-Beyer, 
Lowrance-Wagner, Longest Common Subsequence, and Shortest Common Supersequence algorithms. 

***Please ensure that numpy is installed so that this project can work correctly***

# Installation and Usage
> [!IMPORTANT]
> Not every algorithm uses every method.
> Please refer to implementation table to see which methods each algorithm can perform.

```
pip install goombay
```

All algorithms have a class with customizable parameters and a class instance with default parameters.

Each algorithm is able to perform tasks such as alignment and displaying the underlying matrices, as shown in the implementation table. All algorithms, with the exception of the Hirschberg algorithm, can perform distance, similarity, normalized distance, and normalized similarity calculations.

The methods for the algorithms are:

1. `.distance(seq1, seq2)` - integer value representing the distance between two sequences based on **match score**, **mismatch penalty**, and **gap penalties**.

2. `.similarity(seq1, seq2)` - integer value representing the similarity between two sequences based on **match score**, **mismatch penalty**, and **gap penalties**.

3. `.normalized_distance(seq1, seq2)` - float between `0` and `1`; `0` representing two identical sequences and `1` representing two sequences with no similarities.

4. `.normalized_similarity(seq1, seq2)` - float between `0` and `1`; `1` representing two identical sequences and `0` representing two sequences with no similarities.

5. `.align(seq1, seq2)` - formatted string of the alignment between the provided sequences.

6. `.matrix(seq1, seq2)` - matrix (or matrices) created through the dynamic programming process.

The Hamming distance has two additional methods called `.binary_distance_array` and `.binary_similarity_array` that produce a list of bits denoting which pairwise combinations are a match and which are a mismatch.

# Implementation

**Below is a table of the methods implemented for each algorithm as well as the class (customizable) and instance (default parameters) names.**

| Algorithm                    | Alignment | Matrices | Distance/Similarity/Normalized | Class                         | Instance                      |
| ------------------           | --------- | -------- | ------------------------------ | ----------------------------- | ----------------------------- |
|Needleman-Wunsch              |    [x]    |    [x]   |               [x]              | Needleman_Wunsch              | needleman_wunsch              |
|Gotoh (Global)                |    [x]    |    [x]   |               [x]              | Gotoh                         | gotoh                         |
|Gotoh (Local)                 |    [x]    |    [x]   |               [x]              | Gotoh_Local                   | gotoh_local                   |
|Smith-Waterman                |    [x]    |    [x]   |               [x]              | Smith_Waterman                | smith_waterman                |
|Waterman-Smith-Beyer          |    [x]    |    [x]   |               [x]              | Waterman_Smith_Beyer          | waterman_smith_beyer          |
|Wagner-Fischer                |    [x]    |    [x]   |               [x]              | Wagner_Fischer                | wagner_fischer                |
|Lowrance-Wagner               |    [x]    |    [x]   |               [x]              | Lowrance_Wagner               | lowrance_wagner               |
|Hamming                       |    [x]    |    [ ]   |               [x]              | Hamming                       | hamming                       |
|Hirschberg                    |    [x]    |    [ ]   |               [ ]              | Hirschberg                    | hirschberg                    |
|Jaro                          |    [ ]    |    [x]   |               [x]              | Jaro                          | jaro                          |
|Jaro Winkler                  |    [ ]    |    [x]   |               [x]              | Jaro_Winkler                  | jaro_winkler                  |
|Longest Common Subsequence    |    [x]    |    [x]   |               [x]              | Longest_Common_Subsequence    | longest_common_subsequence    |
|Shortest Common Supersequence |    [x]    |    [x]   |               [x]              | Shortest_Common_Supersequence | shortest_common_supersequence |

## Algorithms Explained
- [Hamming](https://en.wikipedia.org/wiki/Hamming_distance) -
  The Hamming distance is a distance measurement between two sequences of the same length which measures the minimum number of substitutions 
  needed to convert one string into the other.
  When comparing numbers, the hamming distance first converts the numbers into binary and then determines the minimum number of bits that need to be flipped to turn
  one binary sequence into the other.
  The implementation in this project measures sequences of different lengths by comparing the letters of the longer sequence against the blanks of the shorter sequence.
  
- [Wagner-Fischer](https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm) - **Levenshtein distance** -
  The Wagner-Fischer algorithm is a global alignment algorithm that computes the Levenshtein distance between two sequences.
  This algorithm has an invariable gap penalty of 1 and a mismatch (or substitution) cost of 1. Matches are worth 0 therefore they do not affect the score.

- [Lowrance-Wagner](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2819-0) - **Damerau–Levenshtein distance**
  The Lowrance-Wagner algorithm is a global alignment algorithm that computes the Levenshtein distance between two sequences 
  with the addition of adjacent swapping between matching adjacent characters.
  

- [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) -
  The Needleman-Wunsch algorithm is a global alignment algorithm that uses a generalized form of the Levenshtein distance 
  which allows for different weights to be given to matches, mismatches, and gaps.
  - The keyword arguments for this algorithm are `match_score:int = 0`, `mismatch_penalty:int = 1`, and `gap_penalty:int = 2`.

- [Gotoh (Global)](https://helios2.mi.parisdescartes.fr/~lomn/Cours/BI/Material/gap-penalty-gotoh.pdf) -
  The Gotoh algorithm is a global alignment algorithm that is a modification to the Levenshtein distance that uses an affine gap penalty
  (similar to the Waterman-Smith-Beyer algorithm)
  that differentiates between newly created gaps and continuations of gaps.
  This algorithm uses three matrices; ***D*** (optimal score under affine gap penalties), ***P*** (optimal score given that query sequence ends in a gap), and 
  ***Q*** (optimal score given that subject sequence ends in a gap).
  - The keyword arguments for this algorithm are `match_score:int = 0`, `mismatch_penalty:int = 1`, `new_gap_penalty:int = 2`, and `continue_gap_penalty: int = 1`.

- [Gotoh (Local)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh%20(Local)) -
  Similar to the global alignment version of the Gotoh alignment algorithm, the local alignment version also uses three matrices.
  The primary difference is that the optimal alignment score is chosen between applying a penalty for either a mismatch or gap, adding to the total for a match, or zero.
  This allows the cell to be reset to zero if it were to become negative.
  - The keyword arguments for this algorithm are `match_score:int = 1`, `mismatch_penalty:int = 1`, `new_gap_penalty:int = 2`, and `continue_gap_penalty: int = 1`.

- [Smith-Waterman ](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) -
  The Smith-Waterman algorithm is the local alignment equivalent to the Needleman-Wunsch algorithm. Similar to Needleman-Wunsch, it generalizes the Levenshtein distance.
  Similar to the Gotoh local algorithm, it resets any negative cell to zero.
  - The keyword arguments for this algorithm are `match_score:int = 1`, `mismatch_penalty:int = 1`, and `gap_penalty:int = 2`.

- [Waterman-Smith-Beyer](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer) -
  The Waterman-Smith-Beyer algorithm is a global alignment algorithm that is a modification to the Levenshtein distance which uses an arbitrary gap-scoring method.
  The specific implementation used in this package is the affine gap penalty.
  However, a logarithmic or a quadratic gap calculation can also be performed.
  - The keyword arguments for this algorithm are `match_score:int = 0`, `mismatch_penalty:int = 1`, `new_gap_penalty:int = 1`, and `continue_gap_penalty:int = 1`.

- [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm) -
  The Hirschberg algorithm is intended to improve the Needleman-Wunsch algorithm by using recursion to improve space efficiency.
  It uses a method known as divide and conquer to compare the two sequences.
  - The keyword arguments for this algorithm are `match_score: int = 1`, `mismatch_penalty: int = -1`, and `gap_penalty: int = -2`.

- [Jaro & Jaro-Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) -
  The Jaro algorithm is a global alignment algorithm that measures the Jaro distance between two sequences. It produces a number between 0 and 1 that accounts
  for the length of the strings, the number of matching characters, and the number of transpositions. The Jaro algorithm also takes into consideration matches
  that are a certain distance away ((max sequence length/2)-1). Of these matches, transpositions (matches that aren't in the right order) are factored in.

  The Jaro-Winkler algorithm is the same as the Jaro algorithm but also favors sequences that have matching prefix characters (up to four) and adds a scaling factor.
  - The keyword argument for the Jaro-Winkler algorithm is `scaling_factor = 0.1`. The scaling factor should not exceed 0.25 or else it may be possible for the similarity score to be greater than 1.

- [Longest Common Subsequence](https://en.wikipedia.org/wiki/Longest_common_subsequence) -
  The Longest Common Subsequence algorithm generates a distance score by only allowing deletes while not changing the relative order of the characters.
  This will display all of the shared characters between the sequences.

- [Shortest Common Supersequence](https://en.wikipedia.org/wiki/Shortest_common_supersequence) -
  The Shortest Common Supersequence is the shortest combination of the two sequences that contains all the characters within both sequences
  and does not change the relative order of the characters.


# Code Examples

**Hamming Distance**
```python
from goombay import hamming

qs = "AFTG"
ss = "ACTG"

print(hamming.distance(qs, ss))
# 1
print(hamming.similarity(qs, ss))
# 3
print(hamming.binary_distance_array(qs, ss))
# [0,1,0,0]
print(hamming.binary_similarity_array(qs, ss))
# [1,0,1,1]
print(hamming.normalized_distance(qs, ss))
# 0.25
print(hamming.normalized_similarity(qs, ss))
# 0.75
```

**Needleman-Wunsch**
```python
from goombay import needleman_wunsch

print(needleman_wunsch.distance("ACTG","FHYU"))
# 4
print(needleman_wunsch.distance("ACTG","ACTG"))
# 0
print(needleman_wunsch.similarity("ACTG","FHYU"))
# 0
print(needleman_wunsch.similarity("ACTG","ACTG"))
# 4
print(needleman_wunsch.normalized_distance("ACTG","AATG"))
#0.25
print(needleman_wunsch.normalized_similarity("ACTG","AATG"))
#0.75
print(needleman_wunsch.align("BA","ABA"))
#-BA
#ABA
print(needleman_wunsch.matrix("AFTG","ACTG"))
[[0. 2. 4. 6. 8.]
 [2. 0. 2. 4. 6.]
 [4. 2. 1. 3. 5.]
 [6. 4. 3. 1. 3.]
 [8. 6. 5. 3. 1.]]
 ```

# Caveats
> [!CAUTION]
> There are some issues with alignment to be tackled in later releases.

Due to the recursive nature of the Hirschberg algorithm, if a distance score or matrix is needed it is best to use the Needleman-Wunsch algorithm instead.

Note that due to the fact that the Hamming distance does not allow for insertions or deletions, the "aligned sequence" that is returned is just the original sequences in a formatted string. 
This is due to the fact that actually aligning the two sequences using this algorithm would just lead to two lines of the query sequence. 
It should also be noted that the Hamming distance is intended to only be used with sequences of the same length. 
To compensate for strings of differing lengths, my algorithm adds 1 extra point to the distance for every additional letter in the longer sequence since this can be seen as "swapping" the space for a letter or vice versa. However, any distance obtained this way **will not reflect an accurate Hamming distance**.

My Waterman-Smith-Beyer implementation does not always align with that of [Freiburg University](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer), the site I've been using for alignment validation.
Their implementation may have an issue and not mine but I wanted to mention this here and provide the link to my [StackOverflow](https://bioinformatics.stackexchange.com/questions/22683/waterman-smith-beyer-implementation-in-python) question for the sake of posterity.

At the beginning of this project, I thought that the Levenshtein distance was an algorithm, but it is the end result that is being calculated with an approach such as Wagner-Fischer which uses Needleman-Wunsch-esque matrices to calculate the Levenshtein distance.
Thus, the Levenshtein distance implementation has been switched with the Wagner-Fischer algorithm.
Damerau-Levenshtein distance is found using the Lowrance-Wagner algorithm.

Will have to do some experimenting but it appears that the normalized distance/similarity results have undefined behaviour if the match score is not 0.

For the following sequences
```
    qqs = "AGCTCATCAGTCATGCATCCT"
    sss = "CAG"
```
The Gotoh algorithm produces a suboptimal alignment of 
```
AGCTCATCAGTCATGCATCCT
---------------CAG---
```
Correct alignment should be
```
AGCTCATCAGTCATGCATCCT
-------CAG-----------
```

