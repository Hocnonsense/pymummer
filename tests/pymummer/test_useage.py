# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 17:25:29
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-13 21:51:50
 * @FilePath: /pymummer/tests/pymummer/test_useage.py
 * @Description:
"""
# """

from io import StringIO

from pymummer import delta, usage
from tests import Path, temp_output, test_files, test_temp


def test_doc():
    import doctest

    doctest.testmod(usage, raise_on_error=True)


delta_file = test_files / "MarsFilter2.delta"


def test_delta_drop():
    d = usage.Delta_drop(delta_file)


def test_delta_cov():
    d = delta.Delta(delta_file, {})
    with StringIO() as buf:
        s = usage.report_flattern_cov(d.flattern["ref"], d.query.stem).to_csv(buf)
        buf.seek(0)
        assert buf.read() == (
            "Contig,NODE_1564_length_766_cov_111365.326301,NODE_1564_length_766_cov_111365.326301,NODE_1652_length_710_cov_139106.876336,NODE_1652_length_710_cov_139106.876336,NODE_733_length_1545_cov_136201.359060,NODE_733_length_1545_cov_136201.359060\n"
            "Breadth,1,2,0,1,0,1\n"
            "Label,,,,,,\n"
            "A501-plasmid,725.0,41.0,1.0,709.0,48.0,1497.0\n"
        )
    with StringIO() as buf:
        s = usage.report_flattern_cov(d.flattern["query"], d.ref.stem).to_csv(buf)
        buf.seek(0)
        assert buf.read() == (
            "Contig,NZ_CP008888.1,NZ_CP008888.1,NZ_CP008888.1\n"
            "Breadth,0,1,2\n"
            "Label,,,\n"
            "MarsFilter2-sub,843.0,2558.0,228.0\n"
        )


def test_delta_column():
    d = delta.Delta(delta_file, {})
    flattern_align = d.flattern["query"]
    with StringIO() as buf, open(test_files / "compare" / "test_delta_str.tsv") as ref:
        usage.report_flattern_diff(
            flattern_align, min_diff=3, include_unaligned=True, stdout=buf
        )
        buf.seek(0)
        assert buf.read() == ref.read()
    with StringIO() as buf:
        usage.report_flattern_diff(flattern_align, min_diff=3, stdout=buf)
        buf.seek(0)
        assert buf.read() == (
            "#>SeqID\n"
            "#loc\tn_identical\tn_diff\tmin_diff\t*aligns\n"
            ">NZ_CP008888.1\n"
            "3571\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3572\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3573\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3574\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3575\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3576\t0\t1\t1\tNODE_1564_length_766_cov_111365.326301[0:238](+): C\n"
            "3577\t0\t1\t1\tNODE_1564_length_766_cov_111365.326301[0:238](+): T\n"
            "3578\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3579\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3580\t0\t1\t1\tNODE_1564_length_766_cov_111365.326301[0:238](+): A\n"
            "3581\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3582\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3583\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3584\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3585\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3586\t0\t1\t1\tNODE_1564_length_766_cov_111365.326301[0:238](+): G\n"
            "3587\t0\t1\t1\tNODE_1564_length_766_cov_111365.326301[0:238](+): G\n"
            "3588\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
            "3589\t1\t0\t0\tNODE_1564_length_766_cov_111365.326301[0:238](+): |\n"
        )
