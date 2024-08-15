# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 17:23:26
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-15 20:16:02
 * @FilePath: /pymummer/pymummer/usage.py
 * @Description:
"""
# """

from pathlib import Path
from sys import stdout

from .alignment import AlignRegion
from .delta import Delta, DeltaRegion
from .pair import Pair


def Delta_drop(delta_file: Path, cache: dict | bool = True) -> Delta:
    cache_ = ({} if cache else None) if isinstance(cache, bool) else cache
    d = Delta(delta_file, cache_)
    d.drop_alter("ref")
    return d


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
                    i_merged = i_merged.merge(i)
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
