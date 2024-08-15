# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-11 21:02:48
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-15 17:45:47
 * @FilePath: /pymummer/tests/pymummer/test_alignment.py
 * @Description:
"""
# """

from tests import Path, temp_output, test_files, test_temp

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

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
    assert ar2.alignment is not None
    assert "-||+|||-||||" == ar1.alignment == ar3.alignment
    assert "-|+||||-||||" == ar2.alignment[::-1]
    diff = list(ar1.diff or ())
    assert ["-", "|", "|", "g+|", "|", "|", "-", "|", "|", "|", "|"] == diff
    assert "".join(
        [
            bdiff.replace("-", "").replace("+", "").replace("|", bref)
            for bref, bdiff in zip(str(ar1.seq["ref"]), diff)
        ]
    ) == str(ar1.seq["query"])


def test_flattern():
    ac, ar1 = _make_example()
    ar3 = ac.align((SimpleLocation(0, 4, 1), SimpleLocation(0, 5, 1)))
    ac.alignregions.append(ar3)
    flattern = alignment.flatten([ac])
    assert str(ar3.seq["query"]) == "cggta"
    assert flattern == {
        "test1": [
            [(ar1, "-"), (ar3, "-")],  # a -> "", ""
            [(ar1, "|"), (ar3, "|")],  # c -> "c", "c"
            [(ar1, "|"), (ar3, "|")],  # g -> "g", "g"
            [(ar1, "g+|"), (ar3, "g+|+a")],  # t -> "gt", "gta" -> "cggta"
            [(ar1, "|")],  # a
            [(ar1, "|")],  # g
            [(ar1, "-")],  # c -> ""
            [(ar1, "|")],  # t
            [(ar1, "|")],  # g
            [(ar1, "|")],  # a
            [(ar1, "|")],  # g
        ]
    }
    feat = SeqFeature(ar1.loc2["ref"])
    ac.seq2["ref"].id = ac.seqid2["ref"]
    feat2rec = {k: v for k, v in alignment.flatten2feat(feat, flattern, ac.seq2["ref"])}
    assert feat2rec[ar1].seq == ar1.seq["query"]
    assert feat2rec[ar1].annotations["Support"] == [(0, 10)]
    assert str(feat2rec[ar3].seq) == "cggtaagctgag"
    sli: list[tuple[int, int]] = feat2rec[  # pyright: ignore[reportAssignmentType]
        ar3
    ].annotations["Support"]
    assert len(sli) == 1
    assert feat2rec[ar3].seq[slice(*sli[0])] == ar3.seq["query"]
