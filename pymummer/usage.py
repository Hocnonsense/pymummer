# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-12 17:23:26
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-11-17 15:25:14
 * @FilePath: /pymummer/pymummer/usage.py
 * @Description:
"""
# """

from pathlib import Path
from sys import stdout
from typing import Iterable
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from .flatten import flatten2feat, try_get_cds_end
from .alignment import AlignContig2, AlignRegion
from .delta import Delta, DeltaRegion
from .pair import Pair


def Delta_drop(delta_file: Path, cache: dict | bool = True) -> Delta:
    cache_ = ({} if cache else None) if isinstance(cache, bool) else cache
    d = Delta(delta_file, cache_)
    d.drop_alter("ref")
    return d


def report_indel_redund(d: Delta, stdout=stdout):
    write = lambda *x: print(*x, file=stdout)
    this: Pair.T = "ref"
    other = Pair.switch(this)
    ni = 0
    i: list[AlignRegion]
    for g in d.pairs:
        for i in g.dup_dels:  # pyright: ignore[reportAssignmentType]
            if ni:  # new line between each item
                write()
            g.write_mask_regions(i, stdout)
            ni += 1
    return ni


def report_indel_long(
    d: Delta, min_diff=5, skip_aln: set[AlignRegion] | None = None, stdout=stdout
):
    write = lambda *x: print(*x, file=stdout)
    this: Pair.T = "ref"
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
    this: Pair.T = "ref"
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
                # i in g1 should not be redundant -- at least 50? 80? different area
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


def read_prodigal_gff(gff: Path):
    with open(gff) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith(">"):
                break
            if line.startswith("#"):
                continue
            # NODE_1_length_25392_cov_15065.424675 pyrodigal_v3.3.0 CDS 3 1547 179.6 + 0 ID=NODE_1_length_25392_cov_15065.424675_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.410;transl_table=11;conf=99.99;score=179.62;cscore=176.40;sscore=3.22;rscore=0.00;uscore=0.00;tscore=3.22;
            seqid, source, featype, start, end, score, strand, frame, attributes = (
                line.strip().split("\t")
            )
            qualifiers = {
                "source": [source],
                "score": [score],
                "seqid": [seqid],
                "frame": [frame],
            }
            qualifiers.update(
                {
                    k: v.split(",")
                    for k, v in (x.split("=") for x in attributes.strip(";").split(";"))
                }
            )
            start_, stop_ = sorted((int(start), int(end)))
            yield seqid, SeqFeature(
                # Convert from GFF 1-based to standard Python 0-based indexing used by
                # BioPython
                SimpleLocation(start_ - 1, stop_, strand=1 if strand == "+" else -1),
                id=qualifiers["ID"][0],
                type=featype,
                qualifiers=qualifiers,
            )


def report_mut_feat_cds(
    d: Delta,
    feats: Iterable[tuple[str, SeqFeature]],
    template: Pair.T = "ref",
    stdout=stdout,
):
    assert d.seqs
    write = lambda *x: print(*x, file=stdout)
    this = template
    other = Pair.switch(this)
    ni = 0
    nj = 0
    flatten = d.flatten[this]
    for seqid, feat in feats:
        seq2ar_rec = {
            rec.seq: (ar, rec)
            for ar, rec in flatten2feat(feat, flatten, d.seqs[this][seqid])
        }
        if not seq2ar_rec:
            continue
        if len(seq2ar_rec) == 1 and feat.extract(d.seqs[this][seqid].seq) in seq2ar_rec:
            continue
        refcds: SeqRecord = feat.extract(  # pyright: ignore[reportAssignmentType]
            d.seqs[this][seqid]
        )
        has_the_same = refcds.seq in seq2ar_rec
        # print(f"#{seqid} {repr(feat)} {loc} {has_the_same=}")
        if ni + nj:
            write()
        for ar, rec in seq2ar_rec.values():
            assert ar.contig is not None
            if rec.seq == refcds.seq:
                continue
            write(repr(rec), ar, ar.contig)
            ni += 1
            desc = ""
            mutid = f"{feat.id}-{ar.contig.seqid2[template]}-{ar}"
            ar_cds = AlignContig2(
                ("template transcript", "mutation"), (refcds, rec)
            ).align()
            mut_translate = rec.translate()
            refaa = refcds.translate()
            ar_aa = AlignContig2(
                ("template protein", "mutation"), (refaa, mut_translate)
            ).align()
            _feat, new_rec = try_get_cds_end(
                feat, flatten, d.seqs[this][seqid], {ar}, True
            )
            assert ar_cds.contig is not None and ar_aa.contig is not None
            assert rec.seq is not None
            if len(rec.seq) % 3:
                write("Frame shift detected!")
                if "*" in mut_translate.seq[:-1]:
                    desc = "Early terminal!"
                elif "*" not in mut_translate.seq:
                    desc = "Late terminal!"
                if _feat == feat and not rec.seq.endswith("*"):
                    desc += " cannot find a stop codon!"
                ar_cds.contig.write_mask_regions([ar_cds], stdout)
                ar_aa.contig.write_mask_regions([ar_aa], stdout)
                write(f">{mutid} {desc}\n{new_rec.seq}")
            elif len(refcds) == len(rec):
                write("SNP detected!")
                ar_cds.contig.write_mask_regions([ar_cds], stdout)
                if refaa.seq == mut_translate.seq:
                    write("Synonymous mutation (silent)")
                elif len(refaa) == len(mut_translate.seq):
                    write("Non-synonymous mutation (missense)!")
                    if ar_aa.seq and "*" in mut_translate.seq[:-1]:
                        write("Early terminal!")
                    ar_aa.contig.write_mask_regions([ar_aa], stdout)
                    write(f">{mutid} {desc}\n{new_rec.seq}")
                else:
                    write("Non-synonymous mutation (nonsense)!")
                    ar_aa.contig.write_mask_regions([ar_aa], stdout)
                    write(f">{mutid} {desc}\n{new_rec.seq}")
            else:
                write("Indel mutation without frame shift")
                ar_cds.contig.write_mask_regions([ar_cds], stdout)
                ar_aa.contig.write_mask_regions([ar_aa], stdout)
                write(f">{mutid} {desc}\n{new_rec.seq}")


def report_mut_feat_cds2(
    d: Delta,
    feats: Iterable[tuple[str, SeqFeature]],
    template: Pair.T = "ref",
    stdtsv=stdout,
    stdfna=stdout,
):
    assert d.seqs
    writetsv = lambda *x: print(*x, sep="\t", file=stdtsv)
    writetsv("#CHROM", "FEAT", "ALIGN", "TYPE", "HASREF", "HGVS_N", "HGVS_P", "MUTID")
    writefa = lambda i, s: print(f">{i}\n{s}", file=stdfna)
    this = template
    other = Pair.switch(this)
    ni = 0
    nj = 0
    flatten = d.flatten[this]
    for seqid, feat in feats:
        seq2ar_rec = {
            rec.seq: (ar, rec)
            for ar, rec in flatten2feat(feat, flatten, d.seqs[this][seqid])
        }
        if not seq2ar_rec:
            continue
        if len(seq2ar_rec) == 1 and feat.extract(d.seqs[this][seqid].seq) in seq2ar_rec:
            continue
        refcds: SeqRecord = feat.extract(  # pyright: ignore[reportAssignmentType]
            d.seqs[this][seqid]
        )
        has_the_same = refcds.seq in seq2ar_rec
        for ar, rec in seq2ar_rec.values():
            assert ar.contig is not None
            if rec.seq == refcds.seq:
                continue
            ni += 1
            desc = ""
            alnid = f"{ar.contig.seqid2[other]}-{ar}"
            mutid = f"{feat.id}-{alnid}"
            ar_cds = AlignContig2(
                ("template transcript", "mutation"), (refcds, rec)
            ).align()
            mut_translate = rec.translate()
            refaa = refcds.translate()
            ar_aa = AlignContig2(
                ("template protein", "mutation"), (refaa, mut_translate)
            ).align()
            _feat, new_rec = try_get_cds_end(
                feat, flatten, d.seqs[this][seqid], {ar}, True
            )
            assert ar_cds.contig is not None and ar_aa.contig is not None
            assert rec.seq is not None
            mask_aa_hgvs = False
            if len(rec.seq) % 3:
                desc = "indel shift"
                if "*" in mut_translate.seq[:-1]:
                    desc += " early terminal"
                elif "*" not in mut_translate.seq:
                    desc += " late terminal"
                if _feat == feat and "*" not in ar_aa.seq["query"]:
                    desc += " (no stop codon)"
                mask_aa_hgvs = True
            elif len(refcds) == len(rec):
                desc = "snp"
                if refaa.seq == mut_translate.seq:
                    desc += " (silent)"  # Synonymous mutation
                elif len(refaa) == len(mut_translate.seq):
                    desc += " (missense)"  # Non-synonymous mutation
                    if ar_aa.seq and "*" in mut_translate.seq[:-1]:
                        desc += " early terminal"
                else:
                    desc += " (nonsense)"  # Non-synonymous mutation
            else:
                desc = "indel (no shift)"
            writetsv(
                seqid,
                feat.id,
                alnid,
                desc,
                has_the_same,
                " ".join(str(i) for i in ar_cds.hgvs("ref", "g", "@")),
                (
                    "---"
                    if mask_aa_hgvs
                    else " ".join(str(i) for i in ar_aa.hgvs("ref", "p", "@"))
                ),
                mutid,
            )
            writefa(f"{mutid} {desc}", f"{new_rec.seq}")
