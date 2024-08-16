# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-15 18:20:33
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-16 22:04:56
 * @FilePath: /pymummer/pymummer/flatten.py
 * @Description:
"""

from typing import Collection, Iterable, TypeVar, overload
from sys import stdout

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from . import delta

from .pair import Pair
from .alignment import AlignContig2, AlignRegion


TRegion = TypeVar("TRegion", AlignRegion, "delta.DeltaRegion", covariant=True)


# """
# @overload
# def flatten(
#    acs: Iterable["delta.DeltaContig2"],
#    target: Pair.T = "ref",
#    flatten_dict: dict[str, list[list[tuple["delta.DeltaRegion", str]]]] | None = None,
# ) -> dict[str, list[list[tuple["delta.DeltaRegion", str]]]]: ...
# @overload
# def flatten(
#    acs: Iterable[AlignContig2],
#    target: Pair.T = "ref",
#    flatten_dict: dict[str, list[list[tuple[AlignRegion, str]]]] | None = None,
# ) -> dict[str, list[list[tuple[AlignRegion, str]]]]: ...
def flatten(
    acs: Iterable["AlignContig2 | delta.DeltaContig2"],
    target: Pair.T = "ref",
    flatten_dict: dict[str, list[list[tuple[TRegion, str]]]] | None = None,
):
    """
    Warning: MUST make sure seqid is identical in all AlignContig2 !!!
    """
    this = target
    other = Pair.switch(this)
    if flatten_dict is None:
        flatten_dict = {}
    for g in acs:
        if g.seqid2[this] not in flatten_dict:
            flatten_dict[g.seqid2[this]] = [
                []
                for i in range(len(g.seq2[this]))  # pyright: ignore[reportArgumentType]
            ]
        i: TRegion
        for i in g.alignregions:  # type: ignore[assignment]
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
                flatten_dict[i.contig.seqid2[this]][align_base[0]].append(
                    (i, align_base[1])
                )
    return flatten_dict


def flatten2feat(
    feat: SeqFeature,
    flatten_align: dict[str, list[list[tuple[TRegion, str]]]],
    rec: SeqRecord,
    circular: bool | None = True,
    keep_ars: Collection[TRegion] | None = None,
):
    """
    cilcular: True: enable circular, extend the end to the start
              False: disable circular, if reach the end, raise error
              None: don't circular, if reach the end, stop but not raise error
    """
    assert feat.location is not None
    seqid = rec.id
    assert seqid is not None
    refseq = rec[:0]
    ar2rec = {}
    ar2start_end: dict[TRegion, list[int]] = {}
    locat_iter = iter(
        feat.location if feat.location.strand == 1 else reversed(list(feat.location))
    )
    if circular:
        locat_iter = (i % len(rec) for i in locat_iter)
    elif circular is None:
        locat_iter = (i for i in locat_iter if i < len(rec))
    for basei in locat_iter:
        refbase = rec[basei]
        for ar, s in flatten_align[seqid][basei]:
            if keep_ars is not None and ar not in keep_ars:
                continue
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
        seq = Seq(ar2rec[ar].seq)
        supp_t = [tuple(supp[i * 2 : i * 2 + 2]) for i in range(len(supp) // 2)]
        if feat.location.strand == -1:
            seq = seq.reverse_complement()
            supp_t = [(len(seq) - j, len(seq) - i) for i, j in supp_t[::-1]]
        yield ar, SeqRecord(
            seq,
            id=feat.id,
            description=f"{ar.contig.seqid2['query']} {ar}",
            annotations={"Support": supp_t}  # pyright: ignore[reportArgumentType]
            | feat.qualifiers,
        )


def try_get_cds_end(
    feat: SeqFeature,
    flatten_align: dict[str, list[list[tuple[TRegion, str]]]],
    rec: SeqRecord,
    ars: Collection[TRegion] | None = None,
    circular: bool | None = True,
) -> tuple[SeqFeature, SeqRecord]:
    assert feat.location is not None
    this_rec = rec[:0]
    for _ar, this_rec in flatten2feat(feat, flatten_align, rec, circular or None, ars):
        trans: Seq = this_rec.translate().seq
        if "*" in trans[:-1]:
            return feat, this_rec[: (trans.index("*") + 1) * 3]
        if trans.endswith("*"):
            return feat, this_rec
    assert len(this_rec) > 0
    start, end = int(feat.location.start), int(  # pyright: ignore[reportArgumentType]
        feat.location.end  # pyright: ignore[reportArgumentType]
    )
    for i in range(0, len(rec) - len(feat.location), 300):
        if feat.location.strand == 1:
            new_feat = SeqFeature(
                location=SimpleLocation(end + i, end + i + 300, feat.location.strand),
                type=feat.type,
                id=feat.id,
                qualifiers=feat.qualifiers.copy(),
            )
        else:
            new_feat = SeqFeature(
                location=SimpleLocation(
                    start - i - 300, start - i, feat.location.strand
                ),
                type=feat.type,
                id=feat.id,
                qualifiers=feat.qualifiers.copy(),
            )
        print(repr(new_feat))
        processed = False
        for _ar, new_rec in flatten2feat(
            new_feat, flatten_align, rec, circular or None, ars
        ):
            processed = True
            new_rec_ = this_rec + new_rec
            trans = new_rec_.translate().seq
            if "*" in trans:
                return new_feat, new_rec_[: (trans.index("*") + 1) * 3]
        if not processed:
            new_rec = new_feat.extract(rec)
            new_rec_ = this_rec + new_rec
            trans = new_rec_.translate().seq
            if "*" in trans:
                return new_feat, new_rec_[: (trans.index("*") + 1) * 3]
        this_rec += new_rec  # pyright: ignore[reportPossiblyUnboundVariable, reportOperatorIssue]
    return feat, this_rec


def report_flatten_cov(
    flatten_align: dict[str, list[list[tuple[TRegion, str]]]], Label=""
):
    from collections import Counter

    import pandas as pd

    d = {k: Counter(len(i) for i in v) for k, v in flatten_align.items()}
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


def report_flatten_diff(
    flatten_align: dict[str, list[list[tuple[TRegion, str]]]],
    n_window=10,
    min_diff=3,
    include_unaligned=False,
    stdout=stdout,
):
    write = lambda *x: print(*x, sep="\t", file=stdout)
    write("#" + ">SeqID")
    write("#" + "loc", "n_identical", "n_diff", "min_diff", "*aligns")
    this: Pair.T = "ref"
    for s in flatten_align:
        contig = flatten_align[s][0][0][0].contig
        if contig is not None and contig.seqid2["query"] == s:
            this = "query"
        break
    other = Pair.switch(this)
    for s in flatten_align:
        write(f">{s}")
        d10_used = {i: True for i in range(n_window)}
        d10: dict[int, tuple[int, int, int, int, list[str]]] = {}
        last_miss = basei = -1
        for basei, basesc in enumerate(flatten_align[s]):
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
