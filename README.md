![Static Badge](https://img.shields.io/badge/Project_Name-Goombay-blue)
[![PyPI version](https://img.shields.io/pypi/v/goombay.svg)](https://pypi.python.org/pypi/goombay)
[![License](https://img.shields.io/pypi/l/goombay.svg)](LICENSE)
![GitHub branch check runs](https://img.shields.io/github/check-runs/dawnandrew100/goombay/master)
![PyPI - Downloads](https://img.shields.io/pypi/dm/goombay)

# Goombay
This python project contains several sequence alignment algorithms that can also produce scoring matrices for Needleman-Wunsch, Gotoh, Smith-Waterman, Wagner-Fischer, Waterman-Smith-Beyer, Wagner-Fischer, Lowrance-Wagner, Longest Common Subsequence, and Shortest Common Supersequence algorithms. 

***Please ensure that numpy is installed so that this project can work correctly***

# Installation and Usage

```
pip install goombay
```

All algorithms have both claases with customizable parameters and a class instance with default parameters.

Each algorithm is able to perform tasks such as alignment, and displaying the underlying matrices as is shown in the implementation table. All algorithms are able to perform distance, similarity, normalized distance, and normalized similarity calculations with the exception of the hirschberg algorithm.

The methods for the algorithms are:

1. `.distance(seq1, seq2)` - integer value of distance between two sequences based on **match score**, **mismatch penalty**, and **gap penalties**.

2. `.similarity(seq1, seq2)` - integer value of similarity between two sequences based on **match score**, **mismatch penalty**, and **gap penalties**.

3. `.normalized_distance(seq1, seq2)` - float between `0` and `1` with `0` representing two identical sequences and `1` representing two sequences with no similarities.

4. `.normalized_similarity(seq1, seq2)` - float between `0` and `1` with `1` representing two identical sequences and `0` representing two sequences with no similarities.

5. `.align(seq1, seq2)` - displays a formated string of the alignment between the provided sequences.

6. `.matrix(seq1, seq2)` - displays matrix (or matrices) created by sequences.

The Hamming distance has two additional methods called `.binary_distance_array` and `.binary_similarity_array` that produces a list of bits denoting which pairwise combinations are a match and which are a mismatch.

# Implementation

**Below is a table of the methods implemented for each algorithm as well as the class (cutomizable) and instance (default parameteres) names.**

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
[Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)

[Gotoh (Global)](https://helios2.mi.parisdescartes.fr/~lomn/Cours/BI/Material/gap-penalty-gotoh.pdf)

[Gotoh (Local)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh%20(Local))

[Smith-Waterman ](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)

[Waterman-Smith-Beyer](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer)

[Wagner-Fischer](https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm) <- Levenshtein distance

[Lowrance-Wagner](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2819-0) <- Damerau–Levenshtein distance (Levenshtein distance plus adjacent swapping)

[Hamming](https://en.wikipedia.org/wiki/Hamming_distance)

[Hirschberg](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)

[Jaro & Jaro-Winkler](https://en.wikipedia.org/wiki/Jaro%E2%80%93Winkler_distance)

[Longest Common Subsequence](https://en.wikipedia.org/wiki/Longest_common_subsequence)

[Shortest Common Supersequence](https://en.wikipedia.org/wiki/Shortest_common_supersequence)

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

Due to the recursive nature of the Hirschberg algorithm, if a distance score or matrix is needed it is best to use the Needleman-Wunsch algorithm instead.

Note that due to the fact that the Hamming distance does not allow for insertions, or deletions, the "aligned sequence" that is returned is just the original sequences in a formatted string. 
This is due to the fact that actually aligning the two sequences using this algorithm would just lead to two lines of the query sequence. 
It should also be noted that the Hamming distance is intended to only be used with sequences of the same length. 
To compensate for strings of differing lengths, my algorithm adds 1 extra point to the distance for every additional letter in the longer sequence since this can be seen as "swapping" the empty space for a letter or vice versa. However, any distance obtained this way **will not reflect an accurate Hamming distance**.

My Waterman-Smith-Beyer implementation does not always align with that of [Freiburg University](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer), the site I've been using for alignment validation.
It is possible that their implementation has an issue and not mine but I wanted to mention this here and provide the link to my [StackOverflow](https://bioinformatics.stackexchange.com/questions/22683/waterman-smith-beyer-implementation-in-python) question for the sake of posterity.

During the beginning of this project I thought that the Levenshtein distance was an algorithm, but it is the end result that is being calculated with an approach such as Wagner-Fischer which uses Needleman-Wunsch-esque matrices to calculate the Levenshtein distance.
Thusly, the Levenshtein distance implementation has been switched with the Wagner-Fischer algorithm.
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

