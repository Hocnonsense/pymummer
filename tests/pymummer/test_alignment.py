# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-11 21:02:48
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-12 10:19:07
 * @FilePath: /pymummer/tests/pymummer/test_alignment.py
 * @Description:
"""
# """

from tests import Path, temp_output, test_files, test_temp

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import Seq, SimpleLocation
from pymummer import alignment


def test_doc():
    import doctest

    doctest.testmod(alignment, raise_on_error=True)


def _make_example():
    ac = alignment.AlignContig2(
        ("test1", "test2"),
        (SeqRecord(Seq("acgtagctgag")), SeqRecord(Seq("cggtagtgag"))),
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
    assert "-||+|||-||||" == ar1.alignment == ar3.alignment
    assert "-|+||||-||||" == ar2.alignment[::-1]
    assert list(ar1.diff) == ["-", "|", "|", "g+|", "|", "|", "-", "|", "|", "|", "|"]
