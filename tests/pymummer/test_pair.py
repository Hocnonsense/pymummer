# -*- coding: utf-8 -*-
"""
* @Date: 2024-08-11 16:13:55
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2024-08-18 22:08:49
* @FilePath: /pymummer/tests/pymummer/test_pair.py
* @Description:
"""
# """

from tests import Path, temp_output, test_files, test_temp

from pymummer import pair


def test_doc():
    import doctest

    doctest.testmod(pair, raise_on_error=True)


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

    assert list(d.values()) == [3, 5]
    assert list(d.items()) == [("ref", 3), ("query", 5)]

    try:
        d["X"]  # type: ignore[index]
    except KeyError:
        pass
    else:  # pragma: no cover
        assert not "should raise KeyError"


def test_align_edlib():
    if not pair.IMPORT_AVAIL_EDLIB:  # pragma: no cover
        return
    assert pair.align_edlib("ABD", "ABCD") == ([-3], 1)
    # {'editDistance': 1, 'alphabetLength': 4, 'locations': [(0, 3)], 'cigar': '2=1D1='}
    assert pair.align_edlib("ABCD", "ABD") == ([3], 1)
    # {'editDistance': 1, 'alphabetLength': 4, 'locations': [(0, 2)], 'cigar': '2=1I1='}
    assert pair.align_edlib(
        "AAAAABCCCCDDDEEE",
        "AAAAABKKKBCCCEEE",
    ) == ([-8, -1, -1, 4, 1, 1], 7)


def test_hgvs_from_mut():
    if not pair.IMPORT_AVAIL_HGVS:  # pragma: no cover
        return
    h = pair.hgvs_from_mut(1, "A", "T", "genexx.1")
    assert str(h) == "genexx.1:g.2A>T"
    h = pair.hgvs_from_mut(63, "YI", "S", "protxx.1", "p")
    assert str(h) == "protxx.1:p.64_65delinsS"
