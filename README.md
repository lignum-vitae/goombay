[![Static Badge](https://img.shields.io/badge/Project_Name-Goombay-blue)](https://github.com/lignum-vitae/goombay)
[![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Flignum-vitae%2Fgoombay%2Fmaster%2Fpyproject.toml)](https://github.com/lignum-vitae/goombay/blob/master/pyproject.toml)
[![PyPI version](https://img.shields.io/pypi/v/goombay.svg)](https://pypi.python.org/pypi/goombay)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14903445.svg)](https://doi.org/10.5281/zenodo.14903445)
[![License](https://img.shields.io/pypi/l/goombay.svg)](https://github.com/dawnandrew100/goombay/blob/master/LICENSE)
[![GitHub branch check runs](https://img.shields.io/github/check-runs/lignum-vitae/goombay/master)](https://github.com/lignum-vitae/goombay)

# Goombay

This Python project contains several sequence alignment algorithms that can also produce scoring matrices for Needleman-Wunsch, Gotoh, Smith-Waterman, Wagner-Fischer, Waterman-Smith-Beyer,
Lowrance-Wagner, Longest Common Subsequence, and Shortest Common Supersequence algorithms.

***This project depends on NumPy. Please ensure that it is installed.***

# Installation and Usage

> [!IMPORTANT]
> Not every algorithm uses every method.
> Please refer to implementation table to see which methods each algorithm can perform.
#### Run one of the following commands to install the package

```
pip install goombay
py -m pip install goombay
python -m pip install goombay
python3 -m pip install goombay
```
#### To run all tests from command line, run one of the following commands from the root
```
py -m unittest discover tests
python -m unittest discover tests
python3 -m unittest discover tests
pytest tests\
```

---

All algorithms have a class with customizable parameters and a class instance with default parameters.

Each algorithm is able to perform tasks such as alignment and displaying the underlying matrices.
See the [Implementation](#implementation) table for details about which methods can be used with each algorithm.

The methods for the algorithms are:

1. `.distance(seq1, seq2)` - integer value representing the distance between two sequences based on **match score**, **mismatch penalty**, and **gap penalties**.

2. `.similarity(seq1, seq2)` - integer value representing the similarity between two sequences based on **match score**, **mismatch penalty**, and **gap penalties**.

3. `.normalized_distance(seq1, seq2)` - float between `0` and `1`; `0` representing two identical sequences and `1` representing two sequences with no similarities.

4. `.normalized_similarity(seq1, seq2)` - float between `0` and `1`; `1` representing two identical sequences and `0` representing two sequences with no similarities.

5. `.align(seq1, seq2)` - formatted string of the alignment between the provided sequences.

6. `.matrix(seq1, seq2)` - matrix (or matrices) created through the dynamic programming process.

The Hamming distance has two additional methods called `.binary_distance_array` and `.binary_similarity_array` that produce a list of bits denoting which pairwise combinations are a match and which are a mismatch.

---

The `scoring_matrix` keyword argument accepts a substitution matrix from the [Biobase](https://github.com/lignum-vitae/biobase) package.

The following algorithms accept the `scoring_matrix` keyword argument as a parameter:

- NeedlemanWunsch
- WatermanSmithBeyer
- Hirschberg
- FengDoolittle (only applies to above mentioned pairwise algorithms)

# Implementation

**Below is a table of the methods implemented for each algorithm as well as the class (customizable) and instance (default parameters) names.**

| Algorithm                     | Alignment             | Matrices              | Distance/Similarity/Normalized  | Class                       | Instance                      |
| ----------------------------- | --------------------- | --------------------- | ------------------------------- | --------------------------- | ----------------------------- |
|Needleman-Wunsch               |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | Needleman_Wunsch            | needleman_wunsch              |
|Gotoh (Global)                 |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | Gotoh                       | gotoh                         |
|Gotoh (Local)                  |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | GotohLocal                  | gotoh_local                   |
|Smith-Waterman                 |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | SmithWaterman               | smith_waterman                |
|Waterman-Smith-Beyer           |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | WatermanSmithBeyer          | waterman_smith_beyer          |
|Wagner-Fischer                 |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | WagnerFischer               | wagner_fischer                |
|Lowrance-Wagner                |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | LowranceWagner              | lowrance_wagner               |
|Feng-Doolittle                 |<ul><li> [x] </li></ul>|<ul><li> [ ] </li></ul>|     <ul><li> [ ] </li></ul>     | FengDoolittle               | feng_doolittle                |
|Hamming                        |<ul><li> [x] </li></ul>|<ul><li> [ ] </li></ul>|     <ul><li> [x] </li></ul>     | Hamming                     | hamming                       |
|Hirschberg                     |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | Hirschberg                  | hirschberg                    |
|Jaro                           |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | Jaro                        | jaro                          |
|Jaro Winkler                   |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | JaroWinkler                 | jaro_winkler                  |
|Longest Common Subsequence     |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | LongestCommonSubsequence    | longest_common_subsequence    |
|Longest Common Substring       |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | LongestCommonSubstringMSA   | longest_common_substring_msa  |
|Longest Common Substring (MSA) |<ul><li> [x] </li></ul>|<ul><li> [ ] </li></ul>|     <ul><li> [x] </li></ul>     | LongestCommonSubstringMSA   | longest_common_substring_msa  |
|Shortest Common Supersequence  |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | ShortestCommonSupersequence | shortest_common_supersequence |
|LIPNS                          |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | LIPNS                       | lipns                         |
|MLIPNS                         |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | MLIPNS                      | mlipns                        |
|Hamann                         |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | Hamann                      | hamann                        |
|Simple Matching Coefficient    |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | SimpleMatchingCoefficient   | simple_matching_coefficient   |
|Length Ratio                   |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | LengthRatio                 | length_ratio                  |
|Prefix                         |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | Prefix                      | prefix                        |
|Postfix                        |<ul><li> [x] </li></ul>|<ul><li> [x] </li></ul>|     <ul><li> [x] </li></ul>     | Postfix                     | postfix                       |
|Ratcliff Obershelp             |<ul><li> [x] </li></ul>|<ul><li> [ ] </li></ul>|     <ul><li> [x] </li></ul>     | RatcliffObershelp           | ratcliff_obershelp            |

## Algorithms Explained

- [Hamming](https://en.wikipedia.org/wiki/Hamming_distance) -
  The Hamming distance is a distance measurement between two sequences of the same length which measures the minimum number of substitutions 
  needed to convert one string into the other.
  When comparing numbers, the hamming distance first converts the numbers into binary and then determines the minimum number of bits that need to be flipped to turn
  one binary sequence into the other.
  The hamming distance can only be calculated between strings of the same length. Numbers will be padded with 0 to match the longer binary sequence.

- [Wagner-Fischer](https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm) - **Levenshtein distance** -
  The Wagner-Fischer algorithm is a global alignment algorithm that computes the Levenshtein distance between two sequences.
  This algorithm has an invariable gap penalty of 1 and a mismatch (or substitution) cost of 1. Matches are worth 0 therefore they do not affect the score.
  - The keyword argument for the align method for this algorithm is `all_alignments: bool = False`

- [Lowrance-Wagner](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2819-0) - **Damerau–Levenshtein distance**
  The Lowrance-Wagner algorithm is a global alignment algorithm that computes the Levenshtein distance between two sequences 
  with the addition of adjacent swapping between matching adjacent characters.
  Like the Wagner-Fischer algorithm, there's an invariable gap penalty and mismatch penalty of 1. Matches are worth 0 and do not affect the score.
  In addition to these penalties, there's an invariable transposition penalty cost of 1.
  - The keyword argument for the align method for this algorithm is `all_alignments: bool = False`


- [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) -
  The Needleman-Wunsch algorithm is a global alignment algorithm that uses a generalized form of the Levenshtein distance 
  which allows for different weights to be given to matches, mismatches, and gaps.
  - The keyword arguments for the class of this algorithm are `match_score:int = 2`, `mismatch_penalty:int = 1`, and `gap_penalty:int = 2`.
  - The keyword argument for the align method for this algorithm is `all_alignments: bool = False`

- [Gotoh (Global)](https://helios2.mi.parisdescartes.fr/~lomn/Cours/BI/Material/gap-penalty-gotoh.pdf) -
  The Gotoh algorithm is a global alignment algorithm that is a modification to the Levenshtein distance that uses an affine gap penalty
  (similar to the Waterman-Smith-Beyer algorithm)
  that differentiates between newly created gaps and continuations of gaps.
  This algorithm uses three matrices; ***D*** (optimal score under affine gap penalties), ***P*** (optimal score given that query sequence ends in a gap), and 
  ***Q*** (optimal score given that subject sequence ends in a gap).
  - The keyword arguments for the class of this algorithm are `match_score:int = 2`, `mismatch_penalty:int = 1`, `new_gap_penalty:int = 2`, and `continue_gap_penalty: int = 1`.
  - The keyword argument for the align method for this algorithm is `all_alignments: bool = False`

- [Gotoh (Local)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh%20(Local)) -
  Similar to the global alignment version of the Gotoh alignment algorithm, the local alignment version also uses three matrices.
  The primary difference is that the optimal alignment score is chosen between applying a penalty for either a mismatch or gap, adding to the total for a match, or zero.
  This allows the cell to be reset to zero if it were to become negative.
  - The keyword arguments for the class of this algorithm are `match_score:int = 2`, `mismatch_penalty:int = 1`, `new_gap_penalty:int = 3`, and `continue_gap_penalty: int = 2`.

- [Smith-Waterman ](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) -
  The Smith-Waterman algorithm is the local alignment equivalent to the Needleman-Wunsch algorithm. Similar to Needleman-Wunsch, it generalizes the Levenshtein distance.
  Similar to the Gotoh local algorithm, it resets any negative cell to zero.
  - The keyword arguments for this algorithm are `match_score:int = 1`, `mismatch_penalty:int = 1`, and `gap_penalty:int = 2`.

- [Waterman-Smith-Beyer](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer) -
  The Waterman-Smith-Beyer algorithm is a global alignment algorithm that is a modification to the Levenshtein distance which uses an arbitrary gap-scoring method.
  The specific implementation used in this package is the affine gap penalty.
  However, a logarithmic or a quadratic gap calculation can also be performed.
  - The keyword arguments for the class of this algorithm are `match_score:int = 2`, `mismatch_penalty:int = 1`, `new_gap_penalty:int = 4`, and `continue_gap_penalty:int = 1`.
  - The keyword argument for the align method for this algorithm is `all_alignments: bool = False`

- [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm) -
  The Hirschberg algorithm is intended to improve the Needleman-Wunsch algorithm by using recursion to improve space efficiency.
  It uses a method known as divide and conquer to compare the two sequences.
  - The keyword arguments for this algorithm are `match_score: int = 1`, `mismatch_penalty: int = 2`, and `gap_penalty: int = 4`.

- [Feng Doolittle](https://www.cs.auckland.ac.nz/compsci369s1c/lectures/DW-notes/lecture21.pdf) -
  The Feng Doolittle algorithm is a progressive multiple sequence alignment algorithm
  that uses a pairwise implementation such as Needleman Wunsch to determine a distance
  matrix and a clustering algorithm such as the neighbour joining algorithm to determine
  the final alignment.
  - The keyword arguments for the class instantiation are `cluster: str = neighbor_joining` and `pairwise: str = needleman_wunsch`
  - The keyword argument for the align method is `verbose: bool = False`

- [Jaro & Jaro-Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance) -
  The Jaro algorithm is a global alignment algorithm that measures the Jaro distance between two sequences. It produces a number between 0 and 1 that accounts
  for the length of the strings, the number of matching characters, and the number of transpositions. The Jaro algorithm also takes into consideration matches
  that are a certain distance away ((max sequence length/2)-1). Of these matches, transpositions (matches that aren't in the right order) are factored in.

  The Jaro-Winkler algorithm is the same as the Jaro algorithm but also favors sequences that have matching prefix characters (up to four) and adds a scaling factor.
  - The keyword argument for the Jaro-Winkler algorithm is `scaling_factor = 0.1`. The scaling factor should not exceed 0.25 or else it may be possible for the similarity score to be greater than 1.

- [Longest Common Subsequence](https://en.wikipedia.org/wiki/Longest_common_subsequence) -
  The Longest Common Subsequence algorithm generates an alignment by only allowing deletes while not changing the relative order of the characters.
  There may be more than one longest subsequence, in which case, all longest
  subsequences will be shown in a list. Only matches with a length greater than 1 will be shown.

- [Longest Common Substring](https://en.wikipedia.org/wiki/Longest_common_substring) -
  The Longest Common Substring algorithm generates an alignment by finding the longest string that is present in all given sequences.
  Insertions and deletions are not allowed within these strings. There may be
  more than one longest substring, in which case all longest substrings will be shown in a list.
  Only matches with a length greater than 1 will be shown when aligned.
  Matches less than or equal to 1 will return a similarity and distance score of 0.

- [Shortest Common Supersequence](https://en.wikipedia.org/wiki/Shortest_common_supersequence) -
  The Shortest Common Supersequence is the shortest combination of the two sequences that contains all the characters within both sequences
  and does not change the relative order of the characters.

- [Language Independent Product Name Search](http://www.sial.iias.spb.su/files/386-386-1-PB.pdf) -
  LIPNS is an algorithm that determines the similarity between two strings.
  The "is similar" method returns the intended output of the algorithm,
  either True or False (1 or 0 respectively).
  - The keyword argument for LIPNS is `threshold = 0.25`

- [Modfied Language Independent Product Name Search](http://www.sial.iias.spb.su/files/386-386-1-PB.pdf) -
  MLIPNS is a modification to this algorithm that takes into account the number
  of deletions necessary to get from one word to another. A maximum number of deletions
  is allowed before the words are considered to be not similar to each other.
  - The keyword argument for MLIPNS is `max_mismatch = 2`

- [Hamann](https://modbase.compbio.ucsf.edu/pibase/suppinfo/supplementary_info.pdf) -
  The Hamann algorithm is similar to the SMC algorithm, but also takes into account
  mismatches when calculating its score.

- [Simple Matching Coefficient](https://en.wikipedia.org/wiki/Simple_matching_coefficient) -
  The simple matching coefficient is the number of matching attributes divided by
  the total number of attributes. The SMC is also equal to the Hamann similarity
  plus 1 divided by two.

- **Length Ratio** -
  The length ratio algorithm is a simple algorithm that compares the length of the
  shorter sequence to the length of the longer sequence.

- **Prefix** -
  The prefix algorithm measures the number of contiguous matches from the beginning of a string.

- **Postfix** -
  The Postfix algorithm measures the number of contiguous matches from the end of a string.

- [Ratcliff Obershelp](https://en.wikipedia.org/wiki/Gestalt_pattern_matching) -
  The Ratcliff Obershelp algorithm determines the similarity of two or more strings
  by multiplying two by the matching characters and dividing this product by the total
  number of characters in both strings. The matching characters are found by progressive
  iterations of finding the longest common substring.

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
# 0.25
print(needleman_wunsch.normalized_similarity("ACTG","AATG"))
# 0.75
print(needleman_wunsch.align("BA","ABA"))
# -BA
# ABA
print(needleman_wunsch.align("BA","ABA", all_alignments=True))
# ["ACCG\nA-CG", "ACCG\nAC-G"]
print(needleman_wunsch.matrix("AFTG","ACTG"))
[[0. 2. 4. 6. 8.]
 [2. 0. 2. 4. 6.]
 [4. 2. 1. 3. 5.]
 [6. 4. 3. 1. 3.]
 [8. 6. 5. 3. 1.]]
 ```

Needleman-Wunsch with scoring matrix

- Details about scoring matrices can be found at [Biobase](https://github.com/lignum-vitae/biobase)
- Custom scoring matrices also supported

```py
from goombay import NeedlemanWunsch, needleman_wunsch
from biobase.matrix import Blosum


seq1 = "MENSDSLFKLLAEAKGK"
seq2 = "MEQNSDIFKLAQK"

print(needleman_wunsch.align(seq1, seq2))
# ME-NSDSLFKLLAEAKGK
# MEQNSD-IFK-L--A-QK

needleman62 = NeedlemanWunsch(scoring_matrix=Blosum(62))
print(needleman62.align(seq1, seq2))
# ME-NSDSLFKLLAEAKGK
# MEQNSD-IFK-LAQ---K


```

# Contributions

Interested in contributing to Goombay? Please review our [Contribution Guidelines](https://github.com/lignum-vitae/goombay/blob/master/docs/CONTRIBUTING.md) for detailed instructions on how to get involved.

# Caveats

> [!CAUTION]
> Scoring for Gotoh and Water-Smith-Beyer may be incorrect if new gap and continued gap penalties are not adjusted when substitution matrices are in use.
> It is recommended that the the penalties are adjusted so that they are harsher than the worst mismatch possible in your sequence for a give substitution matrix.
> Refer to the `test_SubMatrix.py` file in the `tests` directory for an example.

Note that due to the fact that the Hamming distance does not allow for insertions or deletions, the "aligned sequence" that is returned is just the original sequences in a formatted string.
This is due to the fact that actually aligning the two sequences using this algorithm would just lead to two lines of the query sequence.
It should also be noted that the Hamming distance is intended to only be used with sequences of the same length.

At the beginning of this project, I thought that the Levenshtein distance was an algorithm, but it is the end result that is being calculated with an approach such as Wagner-Fischer which uses Needleman-Wunsch-esque matrices to calculate the Levenshtein distance.
Thus, the Levenshtein distance implementation has been switched to the Wagner-Fischer algorithm.
Damerau-Levenshtein distance is found using the Lowrance-Wagner algorithm.
