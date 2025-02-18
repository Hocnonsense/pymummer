# -*- coding: utf-8 -*-
"""
* @Date: 2024-08-11 21:02:48
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-02-18 21:37:36
 * @FilePath: /pymummer/tests/pymummer/test_alignment.py
* @Description:
"""
# """

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SimpleLocation

from pymummer import alignment, pair

from Bio.SeqRecord import SeqRecord


from tests import Path, temp_output, test_files, test_temp


def test_doc():
    import doctest

    doctest.testmod(alignment, raise_on_error=True)


def _make_example():
    ac = alignment.AlignContig2(
        (SeqRecord(Seq("acgtagctgag"), "test1"), SeqRecord(Seq("cggtagtgag"), "test2"))
    )
    ar1 = ac.align(1)
    ac.alignregions.append(ar1)
    return ac, ar1


def test_align_contig2():
    ac, ar1 = _make_example()
    ar2 = ac.align((-1, -1))
    ar3 = ac.align((SimpleLocation(0, 4, 1), SimpleLocation(0, 5, 1))).merge(
        ac.align((SimpleLocation(6, 11, 1), SimpleLocation(4, 10, 1)))
    )
    # , [1, -3, 4], 3, ac
    assert str(ac) == "AlignContig2(test1[0:11](+), test2[0:10](+))"
    assert str(ar1) == "AlignRegion([0:11](+)~[0:10](+), 11bp-3)"
    assert ac.seq2["ref"].seq == ar1.seq["ref"] == "acgtagctgag"
    assert str(ar1.seq_align["ref"]) == "acg-tagctgag"
    assert ar2.alignment is not None
    assert "-||+|||-||||" == ar1.alignment == ar3.alignment
    assert "-|+||||-||||" == ar2.alignment[::-1]
    diff = list(ar1.diff)
    assert ["-", "|", "|", "g+|", "|", "|", "-", "|", "|", "|", "|"] == diff
    assert "".join(
        [
            bdiff.replace("-", "").replace("+", "").replace("|", bref)
            for bref, bdiff in zip(str(ar1.seq["ref"]), ar1.diff2["ref"])
        ]
    ) == str(ar1.seq["query"])


def test_sub():
    ac = alignment.AlignContig2(
        (
            SeqRecord(Seq("cggtaacgtagctgagacgtagctgaggtgag"), "test1"),
            SeqRecord(Seq("cggtaacgtagctgaggtgag"), "test2"),
        )
    )
    ar1 = ac.align()
    ac.alignregions.append(ar1)
    assert ar1.alignment == "||||||||||||||||--||-|---||-----"
    ar2 = ar1.sub(0, 32, "biopair")
    assert (
        ac.align(align_method="biopair").alignment
        == ar2.alignment
        == "|||||-----------||||||||||||||||"
    )


def test_hgvs():
    if not pair.IMPORT_AVAIL_HGVS:
        return
    ac, ar1 = _make_example()
    assert ["test1:g.1del", "test1:g.3delinsgg", "test1:g.7del"] == [
        str(i) for i in ar1.hgvs("ref", "g")
    ]
