# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-08 20:18:28
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-09-25 10:16:50
 * @FilePath: /pymummer/pymummer/delta.py
 * @Description:
"""
# """


import os
from pathlib import Path
from sys import stdout
from typing import Literal, Sequence, TextIO, overload

from Bio import SeqFeature, SeqIO

from .alignment import AlignContig2, AlignRegion, SeqRecord, SimpleLocation
from . import flatten
from .pair import Pair


class Delta:
    """
    https://github.com/mummer4/mummer/blob/master/docs/nucmer.README
    """

    def __init__(self, delta_file: Path, cache_fa: dict | None = None):
        self.file = delta_file
        with open(delta_file) as fi:
            (self.ref, self.query), self.alignment_type = self._read_file(fi)
            self.pairs = [
                DeltaContig2.from_delta(i, self) for i in self._group_regions(fi)
            ]
        self.seqs: Pair[dict[str, SeqRecord]] = (
            None  # type: ignore[assignment]
            if cache_fa is None
            else Pair(
                {
                    "ref": self._update_cache_fa(str(self.ref), cache_fa),
                    "query": self._update_cache_fa(str(self.query), cache_fa),
                }
            )
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[{self.alignment_type}({self.ref.stem}, {self.query.stem})"

    def __iter__(self):
        for g in self.pairs:
            yield g.alignregion

    @classmethod
    def _update_cache_fa(cls, fa_file: str, cache_fa: dict[str, dict[str, SeqRecord]]):
        if cache_fa and fa_file in cache_fa:
            ref_seq: dict[str, SeqRecord] = cache_fa[fa_file]
        else:
            ref_seq = cache_fa[fa_file] = SeqIO.to_dict(SeqIO.parse(fa_file, "fasta"))
        return ref_seq

    @classmethod
    def _read_file(cls, fi: TextIO):
        """
        The first two lines of the file are
        identical to the .cluster output. The first line lists the two original input
        files separated by a space, and the second line specifies the alignment data
        type, either "NUCMER" or "PROMER".
        """
        ref, query = fi.readline().strip().split()
        alignment_type: Literal["NUCMER", "PROMER"] = fi.readline().strip()  # type: ignore[assignment]
        return ((Path(ref), Path(query)), alignment_type)

    @classmethod
    def _group_regions(cls, fi: TextIO):
        """
        Every grouping of alignment regions have
        a header, just like the cluster's header in the .cluster file. This is a FASTA
        style header and lists the two sequences that produced the following alignments
        after a '>' and separated by a space, after the two sequences are the lengths
        of those sequences in the same order.
        """
        lines: list[str] = []
        for line in fi:
            if line.startswith(">"):
                if lines:
                    yield lines
                lines = []
            lines.append(line)
        if lines:
            yield lines

    def drop_alter(self, aim: Pair.T):
        for g in self.pairs:
            alignments_good: list[DeltaRegion] = []
            alignments_sub: list[DeltaRegion] = []
            last_i = None
            for i in sorted(
                g.alignregions,
                key=lambda x: (
                    x.loc2[aim].start,
                    -x.loc2[aim].end,  # pyright: ignore[reportOperatorIssue]
                    x.n_mismatch,
                ),
            ):
                if last_i is not None:
                    if (
                        last_i.loc2[aim].start
                        <= i.loc2[aim].start  # pyright: ignore[reportOperatorIssue]
                        and i.loc2[aim].end
                        <= last_i.loc2[aim].end  # pyright: ignore[reportOperatorIssue]
                    ):
                        last_reflen = (  # pyright: ignore[reportOperatorIssue]
                            last_i.loc2[aim].end - last_i.loc2[aim].start
                        )
                        reflen = (  # pyright: ignore[reportOperatorIssue]
                            i.loc2[aim].end - i.loc2[aim].start
                        )
                        if (
                            (last_i.n_mismatch + 1) / last_reflen
                            < i.n_mismatch / reflen
                        ) or last_reflen > (reflen + 1):
                            # slightly soft
                            alignments_sub.append(i)
                        else:
                            alignments_good.append(i)
                        continue
                alignments_good.append(i)
                last_i = i
            g.alignregions = alignments_good
            g.alignments_alter = alignments_sub
        return self

    @property
    @Pair
    def flatten(self, this: Pair.T, other: Pair.T):
        flatten_align: dict[str, list[list[tuple[DeltaRegion, str]]]] = {
            s: [[] for i in range(len(q))] for s, q in self.seqs[this].items()
        }
        flatten.flatten(
            self.pairs, this, flatten_align
        )  # pyright: ignore[reportArgumentType]
        return flatten_align

    @overload
    @classmethod
    def run_nucmer(
        cls, ref: Path, query: Path, outprefix: Path, load: Literal[False]
    ) -> Path: ...

    @overload
    @classmethod
    def run_nucmer(
        cls, ref: Path, query: Path, outprefix: Path, load: Literal[True] = True
    ) -> "Delta": ...

    @classmethod
    def run_nucmer(cls, ref: Path, query: Path, outprefix: Path, load=True):
        os.system(f"nucmer {ref} {query} -p {outprefix}")
        outfile = Path((f"{outprefix}.delta"))
        if load:  # pragma: no cover
            return Delta(outfile)
        return outfile


class DeltaContig2(AlignContig2):
    def __init__(
        self,
        seqids: tuple[str, str],
        seqlens: tuple[int, int],
        delta: Delta | None = None,
    ):
        ref_id, query_id = seqids
        ref_len, query_len = seqlens
        self.seqid2 = Pair({"ref": ref_id, "query": query_id})
        self.seqrecordlen = Pair({"ref": ref_len, "query": query_len})
        self.delta: Delta = delta  # type: ignore[assignment]
        self.alignregions: list[  # pyright: ignore[reportIncompatibleVariableOverride]
            DeltaRegion
        ] = []  # type: ignore[assignment]
        self.alignments_alter: list[DeltaRegion] = []

    @property
    def HAS_SEQ(self):
        return self.delta is not None and self.delta.seqs is not None

    @classmethod
    def from_delta(cls, lines: list[str], delta: Delta):
        headline, *alignment = lines
        ref_id, query_id, _ref_len, _query_len = headline[1:].strip().split()
        self = cls((ref_id, query_id), (int(_ref_len), int(_query_len)), delta)
        self.alignregions = [
            DeltaRegion.from_delta(i, self) for i in self._group_regions(alignment)
        ]
        return self

    @property
    @Pair
    def seq2(  # pyright: ignore[reportIncompatibleMethodOverride]
        self, this: "Pair.T", other: "Pair.T"
    ):
        assert self.HAS_SEQ
        return self.delta.seqs[this][self.seqid2[this]]

    @classmethod
    def _get_region_class(cls):
        return DeltaRegion

    @classmethod
    def _group_regions(cls, alignment: list[str]):
        """
           Following this sequence header is the alignment data. Each alignment region
        has a header that describes the start and end coordinates of the alignment in
        each sequence. These coordinates are inclusive and reference the forward strand
        of the current sequence. Thus, if the start coordinate is greater than the end
        coordinate, the alignment is on the reverse strand. The four digits are the
        start and end in the reference sequence respectively and the start and end in
        the query sequence respectively. These coordinates are always measured in DNA
        bases regardless of the alignment data type. The three digits after the starts
        and stops are the number of errors (non-identities), similarity errors (non-
        positive match scores) and non-alpha characters in the sequence (used to count
        stop-codons i promer data).
        """
        prev_0 = 0
        for i in range(len(alignment)):
            if alignment[i].rstrip() == "0":
                yield alignment[prev_0 : i + 1]
                prev_0 = i + 1
        assert prev_0 == len(alignment), f"No 0 found in the {alignment = }"

    @property
    def dup_dels(self):
        alns: dict[tuple, list[DeltaRegion]] = {}
        for i in self.alignregions:
            loc = (
                (i.loc2["ref"].start, i.loc2["ref"].end),
                (i.loc2["query"].start, i.loc2["query"].end),
            )
            alns.setdefault(loc, []).append(i)
        for aln in alns.values():
            if len(aln) > 1:
                yield aln

    def mask_twin_same(
        self, target: Pair.T, aln: Sequence["AlignRegion | DeltaRegion"] | None = None
    ):
        if aln is None:
            for aln in self.dup_dels:  # pyright: ignore[reportAssignmentType]
                break
        assert aln is not None, "No proper aln found"
        return super().mask_twin_same(target, aln)


class DeltaRegion(AlignRegion):
    def __init__(
        self,
        locs: tuple[SimpleLocation, SimpleLocation],
        alignments: list[int],
        reported_mismatch: int | tuple[int, int, int] = (0, 0, 0),
        contig2: "AlignContig2 | None" = None,
    ):
        if isinstance(reported_mismatch, int):
            reported_mismatch = (reported_mismatch, 0, 0)
        super().__init__(locs, alignments, reported_mismatch[0], contig2)
        self.contig = contig2
        self.n_pos, self.n_alpha = reported_mismatch[1:]

    @classmethod
    def from_delta(cls, alignments: list[str], contig2: DeltaContig2 | None = None):
        ref_loc, query_loc, reported_mismatch = cls._parse_header(alignments[0])
        alignments_ = [int(i) for i in alignments[1:]]
        assert alignments_[-1] == 0, f"must end with 0, {alignments = }"
        assert 0 not in alignments_[:-1], f"inner broken {alignments = }"
        return cls((ref_loc, query_loc), alignments_, reported_mismatch, contig2)

    @classmethod
    def _parse_header(cls, header: str):
        """
        Each alignment region
        has a header that describes the start and end coordinates of the alignment in
        each sequence. These coordinates are inclusive and reference the forward strand
        of the current sequence. Thus, if the start coordinate is greater than the end
        coordinate, the alignment is on the reverse strand. The four digits are the
        start and end in the reference sequence respectively and the start and end in
        the query sequence respectively. These coordinates are always measured in DNA
        bases regardless of the alignment data type. The three digits after the starts
        and stops are the number of errors (non-identities), similarity errors (non-
        positive match scores) and non-alpha characters in the sequence (used to count
        stop-codons i promer data).
        """
        ref_start, ref_end, query_start, query_end, *errors = [
            int(i) for i in header.strip().split()
        ]
        assert len(errors) == 3
        if ref_start < ref_end:
            _ref_start, _ref_end = ref_start, ref_end
            _ref_strand = 1
        else:
            _ref_start, _ref_end = ref_end, ref_start
            _ref_strand = -1
        if query_start < query_end:
            _query_start, _query_end = query_start, query_end
            _query_strand = 1
        else:
            _query_start, _query_end = query_end, query_start
            _query_strand = -1
        return (
            SeqFeature.SimpleLocation(_ref_start - 1, _ref_end, strand=_ref_strand),
            SeqFeature.SimpleLocation(
                _query_start - 1, _query_end, strand=_query_strand
            ),
            (errors[0], errors[1], errors[2]),
        )
