# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-15 18:20:33
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-15 20:02:55
 * @FilePath: /pymummer/pymummer/flatten.py
 * @Description:
"""

from typing import Iterable, Mapping, Sequence
from sys import stdout

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from .pair import Pair
from .alignment import AlignContig2, AlignRegion


# """
def flatten(
    acs: Iterable[AlignContig2],
    target: Pair.T = "ref",
    flattern_dict: dict[str, list[list[tuple[AlignRegion, str]]]] | None = None,
):
    """
    Warning: MUST make sure seqid is identical in all AlignContig2 !!!
    """
    this = target
    other = Pair.switch(this)
    if flattern_dict is None:
        flattern_dict = {}
    for g in acs:
        if g.seqid2[this] not in flattern_dict:
            flattern_dict[g.seqid2[this]] = [[] for i in range(len(g.seq2[this]))]
        for i in g.alignregions:
            ref, query = i.seq_align[other], i.seq_align[this]
            assert ref is not None and query is not None
            if i.loc2[this].strand == -1:  # pragma: no cover
                # reverse back
                ref = ref.reverse_complement()
                query = query.reverse_complement()
            for align_base in enumerate(
                i.pair2diff(query, ref),
                i.loc2[this].start,  # pyright: ignore[reportArgumentType]
            ):
                assert i.contig is not None
                flattern_dict[i.contig.seqid2[this]][align_base[0]].append(
                    (i, align_base[1])
                )
    return flattern_dict


def flatten2feat(
    feat: SeqFeature,
    flatten: Mapping[str, Sequence[Sequence[tuple[AlignRegion, str]]]],
    rec: SeqRecord,
):
    assert feat.location is not None
    seqid = rec.id
    assert seqid is not None
    refseq = rec[:0]
    ar2rec = {}
    ar2start_end: dict[AlignRegion, list[int]] = {}
    for basei in feat.location:
        refbase = rec[basei]
        for ar, s in flatten[seqid][basei]:
            if ar not in ar2rec:
                ar2rec[ar] = refseq  # same prefix
                ar2start_end[ar] = [len(refseq), -1]
            elif len(ar2start_end[ar]) % 2:  # pragma: no cover
                ar2start_end[ar][-1] = len(refseq)
                ar2start_end[ar].append(-1)
            ar2rec[ar] += s.replace("-", "").replace("+", "").replace("|", refbase)
            ar2start_end[ar][-1] = -1
        # |   | 1. len(ar2start_end[ar]) % 2
        # |   |       | 2. ar2start_end[ar][-1] == -1
        # | a | True  | True  | new base added (ar2start_end[ar][-1] reset in last iter)
        # | b | True  | False | no new base added (ar2start_end[ar][-1] kept)
        # | c | False | False | no new base added (the last alignment ended, never restart)
        # | c | False | True  | Impossible, what's wrong?
        for ar in ar2start_end:
            if len(ar2start_end[ar]) % 2:
                ar2rec[ar] += refbase
                "disabled as no alignmnet"
            elif ar2start_end[ar][-1] == -1:
                ar2start_end[ar][-1] = len(ar2rec[ar])
            else:
                ar2start_end[ar].append(-1)
                ar2rec[ar] += refbase
        refseq += refbase
    for ar, supp in ar2start_end.items():
        assert ar.contig
        yield ar, SeqRecord(
            Seq(ar2rec[ar].seq),
            id=feat.id,
            description=f"{ar.contig.seqid2['query']} {ar}",
            annotations={  # pyright: ignore[reportArgumentType]
                "Support": [
                    tuple(supp[i * 2 : i * 2 + 2]) for i in range(len(supp) // 2)
                ],
            }
            | feat.qualifiers,
        )


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
