# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 17:25:29
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-13 15:18:40
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


def test_delta_str():
    d = delta.Delta(delta_file, {})
    flattern_align = d.flattern["query"]
    with StringIO() as buf, open(test_files / "compare" / "test_delta_str.tsv") as ref:
        usage.report_flattern_diff(flattern_align, min_diff=3, stdout=buf)
        buf.seek(0)
        assert buf.read() == ref.read()
