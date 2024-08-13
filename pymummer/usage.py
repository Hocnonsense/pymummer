# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 17:23:26
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-13 16:21:10
 * @FilePath: /pymummer/pymummer/usage.py
 * @Description:
"""
# """

from pathlib import Path
from sys import stdout
from .delta import Delta, DeltaRegion
from .pair import Pair


def Delta_drop(delta_file: Path, cache: dict | bool = True) -> Delta:
    cache_ = ({} if cache else None) if isinstance(cache, bool) else cache
    d = Delta(delta_file, cache_)
    d.drop_alter("ref")
    return d


def report_flattern_cov(
    flattern_align: dict[str, list[list[tuple[DeltaRegion, str]]]], Label=""
):
    import pandas as pd
    from collections import Counter

    d = {k: Counter(len(i) for i in v) for k, v in flattern_align.items()}
    return (
        pd.DataFrame(d)
        .pipe(lambda df: df[sorted(df.columns)])
        .T.pipe(lambda df: df[sorted(df.columns)])
        .rename_axis(["Contig"], axis=0)
        .reset_index()
        .melt(id_vars="Contig", var_name="Breadth", value_name="bp")
        .assign(Label=Label)
        .pivot_table(
            values="bp", index="Label", columns=["Contig", "Breadth"], fill_value=0
        )
    )


def report_flattern_diff(
    flattern_align: dict[str, list[list[tuple[DeltaRegion, str]]]],
    n_window=10,
    min_diff=3,
    stdout=stdout,
):
    write = lambda *x: print(*x, sep="\t", file=stdout)
    write("#" + ">SeqID")
    write("#" + "loc", "n_identical", "n_diff", "min_diff", "*aligns")
    this: Pair.T = "ref"
    for s in flattern_align:
        contig = flattern_align[s][0][0][0].contig
        if contig is not None and contig.seqid2["query"] == s:
            this = "query"
        break
    other = Pair.switch(this)
    for s in flattern_align:
        write(f">{s}")
        d10_used = {i: True for i in range(n_window)}
        d10: dict[int, tuple[int, int, int, int, list[str]]] = {}
        last_miss = basei = -1
        for basei, basesc in enumerate(flattern_align[s]):
            d10[basei % n_window] = (
                basei,
                sum("|" in s for ali, s in basesc),
                sum(s != "|" for ali, s in basesc),
                min(
                    list(len(s.replace("|", "").replace("+", "")) for ali, s in basesc)
                    or [0],
                ),
                [
                    (f"{ali.contig.seqid2[other]}{ali.loc2[other]}: {s}")
                    for ali, s in basesc
                    if ali.contig is not None
                ],
            )
            d10_used[basei % n_window] = False
            if (
                sum(d10[i][2] for i in d10 if d10[i] is not None) >= min_diff
                or sum(d10[i][3] for i in d10 if d10[i] is not None) >= min_diff
                or sum(len(d10[i][4]) == 0 for i in d10 if d10[i] is not None)
                >= min_diff
            ):
                for i in range(basei - n_window + 1, basei + 1):
                    shift = i % n_window
                    if d10_used[shift]:
                        continue
                    d10_used[shift] = True
                    if len(d10[shift][4]) == 0:
                        # since here, we have at least min_diff differences
                        if last_miss == -1:
                            last_miss = d10[shift][0]
                        continue
                    if last_miss != -1:
                        write(
                            f"{last_miss}...{d10[shift][0]-1}",
                            0,
                            d10[shift][0] - last_miss,
                        )
                        last_miss = -1
                    write(*d10[shift][:-1], *d10[shift][-1])
        if last_miss != -1:
            write(f"{last_miss}...{basei}", 0, basei - last_miss + 1)
