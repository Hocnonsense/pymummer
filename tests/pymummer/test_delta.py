# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 14:35:01
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-12 16:02:47
 * @FilePath: /pymummer/tests/pymummer/test_delta.py
 * @Description:
"""
# """

from collections import Counter
from io import StringIO

from pymummer import delta
from tests import Path, temp_output, test_files, test_temp


def test_doc():
    import doctest

    doctest.testmod(delta, raise_on_error=True)


delta_file = test_files / "MarsFilter2.delta"


@temp_output
def test_make_delta(test_temp: Path):
    delta.Delta.run_nucmer(
        Path("tests/file/MarsFilter2-sub.fa"),
        Path("tests/file/A501-plasmid.fa"),
        test_temp / "test",
        load=False,
    )
    try:
        delta.Delta.run_nucmer(
            Path("tests/file/MarsFilter2-sub.fa"),
            Path("tests/file/A501-plasmid.fa"),
            test_temp / "test",
        )
    except Exception:  # pragma: no cover
        print("Never mind")
    else:
        print("WuuHuu")


def test_load_delta():
    d = delta.Delta(delta_file)
    assert d.seqs is None
    cache = {}
    d = delta.Delta(delta_file, cache)
    assert len(cache) == 2
    d.pairs[0].seq2["ref"]
    d = delta.Delta(delta_file, {k: {} for k in cache})
    try:
        d.pairs[0].seq2["ref"]
    except KeyError:
        pass
    else:  # pragma: no cover
        assert False


def test_delta_str():
    d = delta.Delta(delta_file, {})
    d.drop_alter("ref")
    with StringIO() as buf:
        print(d, file=buf)
        for i in d:
            print(i.contig, file=buf)
            print(i, file=buf)
        buf.seek(0)
        assert buf.readlines() == [
            "Delta[NUCMER(MarsFilter2-sub, A501-plasmid)\n",
            "DeltaContig2(NODE_1564_length_766_cov_111365.326301[0:238](+), NZ_CP008888.1[3367:3605](-) ..1)\n",
            "DeltaRegion([0:238](+)~[3367:3605](-), 238bp-6)\n",
            "DeltaContig2(NODE_1652_length_710_cov_139106.876336[1:710](+), NZ_CP008888.1[1442:2151](-))\n",
            "DeltaRegion([1:710](+)~[1442:2151](-), 709bp-0)\n",
            "DeltaContig2(NODE_733_length_1545_cov_136201.359060[0:1497](+), NZ_CP008888.1[0:1497](-))\n",
            "DeltaRegion([0:1497](+)~[0:1497](-), 1497bp-0)\n",
        ]


def test_delta_flattern():
    d = delta.Delta(delta_file, {})
    flattern_align = d.flattern["query"]
    assert {1: 2558, 0: 843, 2: 228} == Counter(
        len(i) for i in flattern_align["NZ_CP008888.1"]
    )
    assert {k: Counter(len(i) for i in v) for k, v in d.flattern["ref"].items()} == {
        "NODE_1564_length_766_cov_111365.326301": {1: 725, 2: 41},
        "NODE_1652_length_710_cov_139106.876336": {1: 709, 0: 1},
        "NODE_733_length_1545_cov_136201.359060": {1: 1497, 0: 48},
    }
