# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 17:23:26
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-14 21:45:57
 * @FilePath: /pymummer/pymummer/usage.py
 * @Description:
"""
# """

from pathlib import Path
from sys import stdout
from typing import Mapping, Sequence

from .alignment import AlignRegion
from .delta import Delta, DeltaRegion
from .pair import Pair


def Delta_drop(delta_file: Path, cache: dict | bool = True) -> Delta:
    cache_ = ({} if cache else None) if isinstance(cache, bool) else cache
    d = Delta(delta_file, cache_)
    d.drop_alter("ref")
    return d


def report_flattern_cov(
    flattern_align: Mapping[str, Sequence[Sequence[tuple[AlignRegion, str]]]], Label=""
):
    from collections import Counter

    import pandas as pd

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
    flattern_align: Mapping[str, Sequence[Sequence[tuple[AlignRegion, str]]]],
    n_window=10,
    min_diff=3,
    include_unaligned=False,
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
        d10: dict[int, tuple[int, int, int, int, Sequence[str]]] = {}
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
                or (
                    include_unaligned
                    and sum(len(d10[i][4]) == 0 for i in d10 if d10[i] is not None)
                    >= min_diff
                )
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


def report_indel_redund(d: Delta, stdout=stdout):
    write = lambda *x: print(*x, file=stdout)
    this = "ref"
    other = Pair.switch(this)
    ni = 0
    for g in d.pairs:
        for i in g.dup_dels:
            if ni:  # new line between each item
                write()
            g.write_mask_regions(i, stdout)
            ni += 1
    return ni


def report_indel_long(
    d: Delta, min_diff=5, skip_aln: set[AlignRegion] | None = None, stdout=stdout
):
    write = lambda *x: print(*x, file=stdout)
    this = "ref"
    other = Pair.switch(this)
    ni = 0
    for g in d.pairs:
        for i in g.alignregions:
            if skip_aln and i in skip_aln:
                continue
            if (
                not i.alignment
                or max(len(i) for i in i.alignment.split("|")) < min_diff
            ):
                continue
            if ni:
                write()
            write(g)
            g.write_mask_regions([i], stdout)
            ni += 1
    return ni


def report_indel_looong(
    d: Delta, skip_aln: set[DeltaRegion] | None = None, stdout=stdout
):
    write = lambda *x: print(*x, file=stdout)
    this = "ref"
    other = Pair.switch(this)
    ni = 0
    nj = 0
    for g in d.pairs:
        if len(g.alignregions) <= 1:
            continue
        if list(g.dup_dels) and all(set(g.alignregions) == set(d) for d in g.dup_dels):
            continue
        # sort by the order on the long sequence, first reversed, then pos strand
        alignregion_: list[DeltaRegion] = sorted(  # pyright: ignore[reportCallIssue]
            set(g.alignregions) - (skip_aln or set()),
            key=lambda x: int(
                x.loc2[other].start  # pyright: ignore[reportArgumentType]
                if x.loc2[other].strand == 1
                else -x.loc2[other].end  # pyright: ignore[reportOperatorIssue]
            ),
        )
        regroups: list[list[DeltaRegion]] = []
        prev = None
        for i in alignregion_:
            if prev is None or i.loc2[other].strand != prev.loc2[other].strand:
                regroups.append([i])
            elif (  # pyright: ignore[reportOperatorIssue]
                i.loc2[this].start < prev.loc2[this].start
            ):
                # the short part of the short sequence aligned again
                regroups.append([i])
                if (
                    prev.loc2[other].start < 10  # pyright: ignore[reportOperatorIssue]
                    and i.loc2[this].start < 10  # pyright: ignore[reportOperatorIssue]
                    and prev.loc2[other].end  # pyright: ignore[reportOperatorIssue]
                    < i.loc2[other].start
                    and i.loc2[this].strand == prev.loc2[this].strand
                ):
                    # already prev.ref_start < i.ref_start
                    # this may be the circular situation
                    # DeltaContig2(NODE_4_length_14728_cov_16307.132625[0:472](+), CP046447.1[1588207:1588679](+) ..1)
                    # [ masked ] NODE_4_length_14728_cov_16307.132625 [472:14728](+)
                    #   14256    0 : 0
                    # [bp same ] CP046447.1 [0:14256](+)
                    # [ masked ] NODE_4_length_14728_cov_16307.132625 [0:472](+)
                    #    472     0 : 0
                    # [bp same ] CP046447.1 [1588207:1588679](+)
                    nj += 1
            # elif i.loc2[other].start == prev.loc2[other].start:
            #    pass
            else:
                regroups[-1].append(i)
            prev = i
        if ni:
            write()
        write(g)
        for g1 in regroups:
            if g1 != regroups[0]:
                write()
            if len(g1) == 1:
                g.write_mask_regions(g1, stdout)
            else:
                i_merged = g1[0]
                for i in g1[1:]:
                    i_merged = g1[0].merge(i)
                # i in g1 should not be redundant -- at least 50? 80? differnt area
                if len(i_merged.loc2[this]) < min(len(i.loc2[this]) for i in g1) * 1.5:
                    # please, not simply overlap
                    for i in g1:
                        write(i)
                    continue
                for i in g1:
                    g.write_mask_regions([i], stdout)
                write("try to merge above alignments:")
                g.write_mask_regions(
                    [i_merged], stdout  # pyright: ignore[reportArgumentType]
                )
            ni += 1
    return ni, nj
