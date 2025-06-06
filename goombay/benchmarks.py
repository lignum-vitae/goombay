from random import randint, choice
from string import ascii_lowercase
from timeit import timeit
import goombay as gb

"""
This module displays the run time for the distance metric of each algorithm that has a distance metric
The purpose of this project is not speed but added functionality so this module is only for reference
This module is intended to be run from the command line as follows -> python benchmarks.py
"""

test_strings1 = [
    "".join(choice(ascii_lowercase) for _ in range(randint(50, 100)))
    for _ in range(100)
]
test_strings2 = [
    "".join(choice(ascii_lowercase) for _ in range(randint(50, 100)))
    for _ in range(100)
]


def test_NW():
    [gb.needleman_wunsch.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_SW():
    [gb.smith_waterman.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_WSB():
    [
        gb.waterman_smith_beyer.distance(a, b)
        for a, b in zip(test_strings1, test_strings2)
    ]


def test_WF():
    [gb.wagner_fischer.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_LW():
    [gb.lowrance_wagner.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_H():
    [gb.hamming.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_J():
    [gb.jaro.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_JW():
    [gb.jaro_winkler.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_LCS():
    [
        gb.longest_common_subsequence.distance(a, b)
        for a, b in zip(test_strings1, test_strings2)
    ]


def test_SCS():
    [
        gb.shortest_common_supersequence.distance(a, b)
        for a, b in zip(test_strings1, test_strings2)
    ]


def test_G():
    [gb.gotoh.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def test_GL():
    [gb.gotoh_local.distance(a, b) for a, b in zip(test_strings1, test_strings2)]


def main():
    print(
        "Each of the following tests creates a list comprehension of 100 sequences that are 50 to 100 charachters each"
    )
    print(
        f"Needleman Wunsch test: Time = {timeit('test_NW()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Smith Waterman test: Time = {timeit('test_SW()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Waterman Smith Beyer test: Time = {timeit('test_WSB()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Wagner Fischer test: Time = {timeit('test_WF()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Lowrace Wagner test: Time = {timeit('test_LW()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Hamming test: Time = {timeit('test_H()', globals = globals(), number=1):0.4f}"
    )
    print(f"Jaro test: Time = {timeit('test_J()', globals = globals(), number=1):0.4f}")
    print(
        f"Jaro Winkler test: Time = {timeit('test_JW()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Longest Common Subsequence test: Time = {timeit('test_LCS()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Shortest Common Supersequence test: Time = {timeit('test_SCS()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Gotoh (Global) test: Time = {timeit('test_G()', globals = globals(), number=1):0.4f}"
    )
    print(
        f"Gotoh (Local) test: Time = {timeit('test_GL()', globals = globals(), number=1):0.4f}"
    )


if __name__ == "__main__":
    main()
