# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-11 16:13:55
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-11 17:08:55
 * @FilePath: /pymummer/tests/pymummer/test_pair.py
 * @Description:
"""
# """

from tests import Path, temp_output, test_files, test_temp

from pymummer import pair


def test_doc():
    import doctest

    doctest.testmod(pair)


def test_pair_decorator():
    from pymummer.pair import Pair

    d = Pair({"ref": 3, "query": 5})
    b = d["query"]
    assert b == 5
    b = d.ref
    assert b == 3

    for i in d:
        assert i in ["ref", "query"]
    for i in Pair.ENUM:
        assert i in ["ref", "query"]

    try:
        d["X"]  # type: ignore[index]
    except KeyError:
        pass
    else:
        assert not "should raise KeyError"


def test_align_edlib():
    from pymummer.pair import align_edlib

    assert align_edlib("ABD", "ABCD") == ([-3], 1)
    # {'editDistance': 1, 'alphabetLength': 4, 'locations': [(0, 3)], 'cigar': '2=1D1='}
    assert align_edlib("ABCD", "ABD") == ([3], 1)
    # {'editDistance': 1, 'alphabetLength': 4, 'locations': [(0, 2)], 'cigar': '2=1I1='}
    assert align_edlib(
        "AAAAABCCCCDDDEEE",
        "AAAAABKKKBCCCEEE",
    ) == ([-8, -1, -1, 4, 1, 1], 7)
