# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-15 18:24:56
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-02-13 10:49:31
 * @FilePath: /pymummer/tests/pymummer/test_flatten.py
 * @Description:
"""
# """

from io import StringIO

from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import pytest

from pymummer import alignment, delta, flatten, usage
from tests import Path, temp_output, test_files, test_temp


def test_doc():
    import doctest

    doctest.testmod(flatten, raise_on_error=True)


def _make_example():
    ac = alignment.AlignContig2(
        (SeqRecord(Seq("acgtagctgag")), SeqRecord(Seq("cggtagtgag"))),
        ("test1", "test2"),
    )
    ar1 = ac.align(1)
    ac.alignregions.append(ar1)
    return ac, ar1


def test_flatten():
    ac, ar1 = _make_example()
    ar3 = ac.align((SimpleLocation(0, 4, 1), SimpleLocation(0, 5, 1)))
    ac.alignregions.append(ar3)
    flatten_align = flatten.flatten([ac])
    assert str(ar3.seq["query"]) == "cggta"
    assert flatten_align == {
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
    feat2rec = {
        k: v for k, v in flatten.flatten2feat(feat, flatten_align, ac.seq2["ref"])
    }
    assert feat2rec[ar1].seq == ar1.seq["query"]
    assert feat2rec[ar1].annotations["Support"] == [(0, 10)]
    assert str(feat2rec[ar3].seq) == "cggtaagctgag"
    sli: list[tuple[int, int]] = feat2rec[  # pyright: ignore[reportAssignmentType]
        ar3
    ].annotations["Support"]
    assert len(sli) == 1
    assert feat2rec[ar3].seq[slice(*sli[0])] == ar3.seq["query"]


delta_file = test_files / "MarsFilter2.delta"


def test_delta_drop():
    d = usage.Delta_drop(delta_file)


def test_delta_cov():
    d = delta.Delta(delta_file, {})
    with StringIO() as buf:
        s = flatten.report_flatten_cov(d.flatten["ref"], d.query.stem).to_csv(buf)
        buf.seek(0)
        assert buf.read() == (
            "Contig,NODE_1564_length_766_cov_111365.326301,NODE_1564_length_766_cov_111365.326301,NODE_1652_length_710_cov_139106.876336,NODE_1652_length_710_cov_139106.876336,NODE_733_length_1545_cov_136201.359060,NODE_733_length_1545_cov_136201.359060\n"
            "Breadth,1,2,0,1,0,1\n"
            "Label,,,,,,\n"
            "A501-plasmid,725.0,41.0,1.0,709.0,48.0,1497.0\n"
        )
    with StringIO() as buf:
        s = flatten.report_flatten_cov(d.flatten["query"], d.ref.stem).to_csv(buf)
        buf.seek(0)
        assert buf.read() == (
            "Contig,NZ_CP008888.1,NZ_CP008888.1,NZ_CP008888.1\n"
            "Breadth,0,1,2\n"
            "Label,,,\n"
            "MarsFilter2-sub,843.0,2558.0,228.0\n"
        )


def test_delta_column():
    d = delta.Delta(delta_file, {})
    flatten_align = d.flatten["query"]
    with StringIO() as buf, open(test_files / "compare" / "test_delta_str.tsv") as ref:
        flatten.report_flatten_diff(
            flatten_align, min_diff=3, include_unaligned=True, stdout=buf
        )
        buf.seek(0)
        assert buf.read() == ref.read()
    with StringIO() as buf:
        flatten.report_flatten_diff(flatten_align, min_diff=3, stdout=buf)
        buf.seek(0)
        assert buf.read() == (
            "#>SeqID\n"
            "#loc\tn_identical\tn_diff\tmin_diff\t*aligns\n"
            ">NZ_CP008888.1\n"
            "3381\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3382\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3383\t1\t1\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): C\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3384\t1\t1\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): C\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3385\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3386\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3387\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3388\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3389\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3390\t1\t1\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): T\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3391\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3392\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3393\t1\t1\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): A\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3394\t1\t1\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): G\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3395\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3396\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3397\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3398\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
            "3399\t2\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\tNODE_1564_length_766_cov_111365.326301[197:766](+): |\n"
        )


def test_try_get_cds_end():
    d = delta.Delta(delta_file, {})
    assert d.seqs
    flatten_align = d.flatten["query"]
    feat = SeqFeature(SimpleLocation(80, 308, -1))
    with pytest.warns(BiopythonWarning, match="Partial codon, "):
        rep_seq = feat.extract(d.seqs["query"]["NZ_CP008888.1"].seq)
        assert (
            rep_seq
            == flatten.try_get_cds_end(
                feat, flatten_align, d.seqs["query"]["NZ_CP008888.1"]
            )[1].seq
        )
        assert (
            rep_seq
            == flatten.try_get_cds_end(
                SeqFeature(SimpleLocation(90, 308, -1)),
                flatten_align,
                d.seqs["query"]["NZ_CP008888.1"],
            )[1].seq
        )
        assert (
            rep_seq
            == flatten.try_get_cds_end(
                SeqFeature(SimpleLocation(40, 308, -1)),
                flatten_align,
                d.seqs["query"]["NZ_CP008888.1"],
            )[1].seq
        )
        feat = SeqFeature(SimpleLocation(80, 185, 1))
        rep_seq = feat.extract(d.seqs["query"]["NZ_CP008888.1"].seq)
        assert (
            rep_seq
            == flatten.try_get_cds_end(
                SeqFeature(SimpleLocation(80, 180, 1)),
                flatten_align,
                d.seqs["query"]["NZ_CP008888.1"],
            )[1].seq
        )
