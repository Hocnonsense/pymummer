# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 17:25:29
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-12 20:55:59
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


def test_delta_str():
    d = delta.Delta(delta_file, {})
    flattern_align = d.flattern["query"]
    with StringIO() as buf, open(test_files / "compare" / "test_delta_str.tsv") as ref:
        usage.report_flattern_diff(flattern_align, min_diff=3, stdout=buf)
        buf.seek(0)
        assert buf.read() == ref.read()
