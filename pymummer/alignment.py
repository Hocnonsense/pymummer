# -*- coding: utf-8 -*-
"""
* @Date: 2024-08-11 17:37:59
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-02-13 10:50:27
* @FilePath: /pymummer/pymummer/alignment.py
* @Description:
"""
# """

import os
from typing import Iterable, Literal, Sequence, overload
from sys import stdout

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .pair import Pair, align_edlib
from pymummer import pair


class AlignContig2:
    @overload
    def __init__(self, seqs: tuple[SeqRecord, SeqRecord]): ...
    @overload
    def __init__(self, seqs: tuple[SeqRecord, SeqRecord], seqids: tuple[str, str]): ...
    @overload
    def __init__(self, seqs: tuple[SeqRecord, SeqRecord], *, refid: str): ...
    @overload
    def __init__(self, seqs: tuple[SeqRecord, SeqRecord], *, queryid: str): ...
    def __init__(
        self,
        seqs: tuple[SeqRecord, SeqRecord],
        seqids: tuple[str, str] | None = None,
        *,
        refid: str | None = None,
        queryid: str | None = None,
    ):
        ref_seq, query_seq = seqs
        ref_id, query_id = ref_seq.id, query_seq.id
        if seqids is not None:
            ref_id, query_id = seqids
        if refid is not None:
            ref_id = refid
        if queryid is not None:
            query_id = queryid
        self.seqid2 = Pair({"ref": ref_id, "query": query_id})
        self.__seq2 = Pair({"ref": ref_seq, "query": query_seq})
        self.alignregions: list[AlignRegion] = []

    @property
    def seq2(self):
        return self.__seq2

    @property
    def alignregion(self):
        return self.alignregions[0]

    def __repr__(self) -> str:
        ref_loc = self.alignregion.loc2["ref"] if self.alignregions else "*"
        query_loc = self.alignregion.loc2["query"] if self.alignregions else "*"
        other_aligns = len(self.alignregions) - 1
        _others = f" ..{other_aligns}" if other_aligns > 0 else ""
        return (
            f"{self.__class__.__name__}"
            f"({self.seqid2['ref']}{ref_loc},"
            f" {self.seqid2['query']}{query_loc}{_others})"
        )

    @overload
    def align(
        self,
        query_strand: Literal[1, -1] = 1,
        /,
        align_method: Literal["edlib"] = "edlib",
    ) -> "AlignRegion": ...
    @overload
    def align(
        self,
        strand: tuple[Literal[1, -1], Literal[1, -1]],
        /,
        align_method: Literal["edlib"] = "edlib",
    ) -> "AlignRegion": ...
    @overload
    def align(
        self,
        loc2: tuple[SimpleLocation, SimpleLocation],
        /,
        align_method: Literal["edlib"] = "edlib",
    ) -> "AlignRegion": ...
    def align(
        self,
        loc2: (
            Literal[1, -1]
            | tuple[Literal[1, -1], Literal[1, -1]]
            | tuple[SimpleLocation, SimpleLocation]
        ) = 1,
        /,
        align_method="edlib",
    ):
        # check locs
        if isinstance(loc2, int):
            loc_ref = SimpleLocation(0, len(self.seq2["ref"]), 1)
            loc_query = SimpleLocation(0, len(self.seq2["query"]), loc2)
        elif isinstance(loc2, tuple):
            _ref, _query = loc2
            if isinstance(_ref, int):
                loc_ref = SimpleLocation(0, len(self.seq2["ref"]), _ref)
            else:
                loc_ref = _ref
            if isinstance(_query, int):
                loc_query = SimpleLocation(0, len(self.seq2["query"]), _query)
            else:
                loc_query = _query
        else:
            raise NotImplementedError(f"not implied {loc2 = }")  # pragma: no cover
        # check align_method and align
        if align_method == "edlib":
            ref_seq = loc_ref.extract(self.seq2["ref"].seq)
            query_seq = loc_query.extract(self.seq2["query"].seq)
            alnm, mism = align_edlib(ref_seq, query_seq)
        else:
            raise NotImplementedError(
                f"not implied {align_method = }"
            )  # pragma: no cover
        # return AlignRegion
        return self.region_class((loc_ref, loc_query), alnm, mism, self)

    @classmethod
    def _get_region_class(cls):
        return AlignRegion

    @property
    def region_class(self):
        return self._get_region_class()

    def mask_twin_same(self, target: Pair.T, aln: Sequence["AlignRegion"]):
        this = target
        other = Pair.switch(this)
        _alignment, *_alignments = [i.alignment2[this] for i in aln]
        assert _alignment is not None
        _alignment1 = _alignment.replace("+", "*").replace("-", "/")
        for _alignment2 in _alignments:
            assert _alignment2 is not None
            _alignment1 = (
                aln[0]
                .pair2align(
                    _alignment1.replace("|", "!"), _alignment2.replace("|", "!")
                )
                .replace("!", "+")
                .replace("+", "*")
                .replace("-", "/")
            )
        last_char = None
        start_i = 0
        align_regions: list[tuple[int, int]] = []
        aln_i = (0, "")
        for aln_i in enumerate(_alignment1):
            if aln_i[1] == "|":
                if last_char != "|":
                    start_i = aln_i[0]
            elif last_char == "|":
                if start_i != 0:
                    start_i += 5
                if aln_i[0] - start_i > 15:
                    align_regions.append((start_i, aln_i[0] - 5))
            last_char = aln_i[1]
        if last_char == "|":
            if start_i != 0:
                start_i += 5
            if aln_i[0] - start_i > 10:
                align_regions.append((start_i, aln_i[0] + 1))
        for i in aln:
            other_align, alignment, this_align = (
                i.seq_align[other],
                i.alignment,
                i.seq_align[this],
            )
            assert other_align is not None
            assert alignment is not None
            assert this_align is not None
            for start, end in reversed(align_regions):
                align_mask_str = f" {end-start:^8} "
                if len(align_mask_str) == 8:
                    ref_mask_str = f"[ masked ]"
                    query_mask_str = f"[bp same ]"
                else:
                    ref_mask_str = "[{k:^{v}}]".format(
                        k="masked", v=len(align_mask_str) - 2
                    )
                    query_mask_str = "[{k:^{v}}]".format(
                        k="bp same", v=len(align_mask_str) - 2
                    )
                alignment = alignment[:start] + align_mask_str + alignment[end:]
                this_align = this_align[:start] + ref_mask_str + this_align[end:]
                other_align = other_align[:start] + query_mask_str + other_align[end:]
            yield i, this_align, alignment, other_align

    def write_mask_regions(self, alignregions: Sequence["AlignRegion"], stdout=stdout):
        write = lambda *x: print(*x, file=stdout)
        this: Pair.T = "ref"
        other = Pair.switch(this)
        _other_align = ""
        _i = None
        for _i, _q, _a, _r in self.mask_twin_same(this, alignregions):
            if not _other_align:
                write(_q, self.seqid2[this], _i.loc2[this])
            write(
                _a,
                len(_i.feat2["ref"].qualifiers["ins"]),
                ":",
                -len(_i.feat2["query"].qualifiers["ins"]),
            )
            _other_align = _r
        assert _i is not None, f"empty {alignregions = }"
        write(_other_align, self.seqid2[other], _i.loc2[other])


class AlignRegion:
    def __init__(
        self,
        locs: tuple[SimpleLocation, SimpleLocation],
        alignments: Iterable[int],
        reported_mismatch=-1,
        contig2: "AlignContig2 | None" = None,
    ):
        self.contig = contig2
        ref_loc, query_loc = locs
        self.loc2 = Pair({"ref": ref_loc, "query": query_loc})
        refid = self.contig.seqid2["ref"] if self.contig else "<unknown id>"
        queryid = self.contig.seqid2["query"] if self.contig else "<unknown id>"
        ref_del, query_del = self._parse_alignment(alignments)
        self.feat2: dict[Pair.T, SeqFeature] = {
            "ref": SeqFeature(
                location=ref_loc,
                type="alignment",
                id=refid,
                qualifiers={"ins": query_del},
            ),
            "query": SeqFeature(
                location=query_loc,
                type="alignment",
                id=queryid,
                qualifiers={"ins": ref_del},
            ),
        }
        self.n_mismatch = reported_mismatch

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}"
            f"({self.loc2['ref']}"
            f"~{self.loc2['query']},"
            f" {self.region_len}bp-{self.n_mismatch})"
        )

    @property
    def HAS_SEQ(self):
        return self.contig is not None

    @property
    def region_len(self):
        return max(len(i) for i in self.loc2.values())

    @property
    @Pair
    def seq(self, this: "Pair.T", other: "Pair.T") -> Seq:
        assert self.HAS_SEQ and self.contig is not None
        seq = self.contig.seq2[this]
        return self.loc2[this].extract(seq.seq)  # pyright: ignore[reportReturnType]

    @classmethod
    def _align_indel(cls, seq: Seq, inserts: Iterable[int]) -> Seq:
        """
        >>> AlignRegion._align_indel("acgtagctgag", [4])
        'acg-tagctgag'
        >>> AlignRegion._align_indel("ABD", [3])
        'AB-D'
        """
        prev = 0
        parts = seq[:0]
        for i in enumerate(sorted(inserts), 1):
            this = i[1] - i[0]
            parts += seq[prev:this] + "-"
            prev = this
        if prev < len(seq):
            parts += seq[prev:]
        return parts

    @property
    @Pair
    def seq_align(self, this: "Pair.T", other: "Pair.T"):
        seq = self.seq[this]
        return self._align_indel(seq, self.feat2[this].qualifiers["ins"])

    @classmethod
    def pair2align(cls, ref, query):
        ref, query = str(ref), str(query)
        assert len(ref) == len(query), f"{ref = } != {query = }"
        parts = []
        prev = 0
        len_seq1 = len(ref)
        while prev < len_seq1:
            cp = os.path.commonprefix([ref[prev:], query[prev:]])
            parts.append("|" * len(cp))
            prev += len(cp)
            if len_seq1 <= prev:
                break
            if ref[prev] == "-":
                parts.append("+")
            else:
                parts.append(query[prev])
            prev += 1
        return "".join(parts)

    @property
    @Pair
    def alignment2(self, this: "Pair.T", other: "Pair.T"):
        ref, query = self.seq_align[this], self.seq_align[other]
        return self.pair2align(ref, query)

    @property
    def alignment(self):
        return self.alignment2["ref"]

    @classmethod
    def pair2diff(cls, ref, query):
        """
        >>> # AB-D
        >>> # ||+|
        >>> # ABCD
        >>> list(AlignRegion.pair2diff("AB-D", "ABCD"))
        ['|', '|', 'C+|']
        >>> # AB-DFGH
        >>> # ||+|E|-
        >>> # ABCDEG-
        >>> list(zip("ABDFGH", AlignRegion.pair2diff("AB-DFGH", "ABCDEG-")))
        [('A', '|'), ('B', '|'), ('D', 'C+|'), ('F', 'E'), ('G', '|'), ('H', '-')]
        >>> list(AlignRegion.pair2diff("AB-", "ABC"))
        ['|', '|+C']
        """
        ref, query = str(ref), str(query)
        assert len(ref.replace("-", "")) > 0, f"empty {ref = }"
        if ref.endswith("-"):
            ref_ = ref.rstrip("-")
            y = cls.pair2diff(ref_, query[: len(ref_)])
            last_y1 = next(y)
            for yi in y:
                yield str(last_y1)
                last_y1 = yi
            yield str(last_y1 + "+" + query[len(ref_) :])
            return
        align = cls.pair2align(ref, query)
        ins = ""
        # print(ref, align, query, sep="\n")
        for base_ref, base_align, base_query in zip(ref, align, query):
            if base_ref == "-":
                ins += base_query
            else:
                if ins:
                    ins += "+"
                if base_align == "|":
                    yield ins + "|"
                else:
                    yield ins + base_query
                ins = ""
        assert not ins, "unexpected base left "

    @property
    @Pair
    def diff2(self, this: "Pair.T", other: "Pair.T"):
        ref, query = self.seq_align[this], self.seq_align[other]
        return self.pair2diff(ref, query)

    @property
    def diff(self):
        return self.diff2["ref"]

    @classmethod
    def diff2vcf(cls, ref, query):
        """
        Unluckly, indel may not be identified by this method.

        >>> # AB-D
        >>> # ||+|
        >>> # ABCD
        >>> list(AlignRegion.diff2vcf("AB-D", "ABCD"))
        [(2, '', 'C')]
        >>> # AB-DFGH
        >>> # ||+|E|-
        >>> # ABCDEG-
        >>> list(AlignRegion.diff2vcf("AB-DFGH", "ABCDEG-"))
        [(2, '', 'C'), (4, 'F', 'E'), (6, 'H', '')]
        >>> list(AlignRegion.diff2vcf("AB-", "ABC"))
        [(2, '', 'C')]
        """
        ref, query = str(ref), str(query)
        assert len(ref.replace("-", "")) > 0, f"empty {ref = }"
        align = cls.pair2align(ref, query)
        _ref = _alt = ""
        index, _i = 0, -1
        for base_ref, base_align, base_query in zip(ref, align, query):
            if base_ref != "-":
                index += 1
            if base_align == "|":
                if _ref or _alt:
                    yield _i, _ref, _alt
                    _ref = _alt = ""
                    _i = -1
            else:
                if _i == -1:
                    _i = index
                if base_ref != "-":
                    _ref += base_ref
                if base_query != "-":
                    _alt += base_query
        if _ref or _alt:
            yield _i, _ref, _alt

    def hgvs(self, this: "Pair.T", seqtype="g", seqid: str | None = None):
        assert pair.IMPORT_AVAIL_HGVS
        if seqid is None:
            assert self.contig
            seqid = self.contig.seqid2[this]
        other = Pair.switch(this)
        ref, query = self.seq_align[this], self.seq_align[other]
        return [
            pair.hgvs_from_mut(s, r, l, seqid, seqtype=seqtype)
            for s, r, l in self.diff2vcf(ref, query)
        ]

    @classmethod
    def mask_muts(cls, alignment: str, n_wing=5, min_mask=10):
        """
        >>> AlignRegion.mask_muts("||||||+|||-|||||", 1, 3)
        [(0, 5), (12, 16)]
        """
        last_char = ""
        start_i = 0
        align_regions: list[tuple[int, int]] = []
        i = 0
        for i, c in enumerate(alignment):
            if c == "|":
                if last_char != "|":
                    # start of a maskable region
                    start_i = i
            elif last_char == "|":
                # end of a maskable region
                if start_i != 0:
                    # if region is at the start,
                    # we don't care about the start sequence
                    # otherwise, don't mask the left n_wing
                    start_i += n_wing
                # don't mask the right n_wing,
                # the inner part is masked, with at least min_mask
                if i - start_i >= min_mask + n_wing:
                    align_regions.append((start_i, i - n_wing))
            last_char = c
        assert i, "Empty alignment!"
        if last_char == "|":
            # end of the alignment
            # also, we don't care about the end sequence
            if start_i != 0:
                start_i += n_wing
            if i - start_i >= min_mask:
                align_regions.append((start_i, i + 1))
        return align_regions

    @classmethod
    def render_mask(
        cls,
        muts: list[tuple[int, int]],
        seq: str,
        label_index=1,
        f_labels: list[str] | tuple[str, str, str] = (
            "[ masked ]",
            " {len:^8} ",
            "[bp same ]",
        ),
    ):
        """
        >>> l = "acgtatctag-tgaga", "||||||+|||-|||||", "acgtat-tagctgaga"
        >>> print(*(AlignRegion.render_mask([(0, 5), (12, 16)], s, i) + " |" for i, s in enumerate(l)), sep="\\n")
        [ masked ]tctag-t[ masked ] |
            5     |+|||-|    4      |
        [bp same ]t-tagct[bp same ] |
        >>> fl = ("[masked]", " {len} ", "[bp same]")
        >>> print(*(AlignRegion.render_mask([(0, 5), (12, 16)], s, i, fl) + " |" for i, s in enumerate(l)), sep="\\n")
        [masked ]tctag-t[masked ] |
            5    |+|||-|    4     |
        [bp same]t-tagct[bp same] |
        >>> fl = ("[masked ", "[ {len:^8} ]", " bp same]")
        >>> print(*(AlignRegion.render_mask([(0, 5), (12, 16)], s, i, fl) + " |" for i, s in enumerate(l)), sep="\\n")
        [masked     tctag-t[masked      |
        [    5     ]|+|||-|[    4     ] |
            bp same]t-tagct    bp same] |
        """
        for start, end in reversed(muts):
            labels = [i.format(len=end - start, v=0) for i in f_labels]
            common_len = max(len(i) for i in labels)
            label = labels[label_index]
            if len(label) < common_len:
                if label.startswith(" "):
                    if label.endswith(" "):
                        label = "{k:^{v}}".format(k=label.strip(), v=common_len)
                    else:
                        label = "{k:>{v}}".format(k=label.strip(), v=common_len)
                elif label.endswith(" "):
                    label = "{k:<{v}}".format(k=label.strip(), v=common_len)
                else:
                    label = "{k0}{k1:^{v}}{k2}".format(
                        k0=label[0], k1=label[1:-1], k2=label[-1], v=common_len - 2
                    )
            seq = seq[:start] + label + seq[end:]
        return seq

    def merge(self, other: "AlignRegion", align_method: Literal["edlib"] = "edlib"):
        assert (
            self.contig == other.contig
        ), f"Not of the same group, this: {self.contig}, other: {other.contig}"
        assert all(
            self.loc2[i].strand == other.loc2[i].strand for i in Pair.ENUM  # type: ignore[index]
        ), f"Strand not the same, this: {self}, other: {other}, contig: {self.contig}"
        assert self.contig is not None

        loc_ref = SimpleLocation(
            min(
                self.loc2["ref"].start,  # pyright: ignore[reportArgumentType]
                other.loc2["ref"].start,  # pyright: ignore[reportArgumentType]
            ),
            max(
                self.loc2["ref"].end,  # pyright: ignore[reportArgumentType]
                other.loc2["ref"].end,  # pyright: ignore[reportArgumentType]
            ),
            self.loc2["ref"].strand,
        )
        loc_query = SimpleLocation(
            min(
                self.loc2["query"].start,  # pyright: ignore[reportArgumentType]
                other.loc2["query"].start,  # pyright: ignore[reportArgumentType]
            ),
            max(
                self.loc2["query"].end,  # pyright: ignore[reportArgumentType]
                other.loc2["query"].end,  # pyright: ignore[reportArgumentType]
            ),
            self.loc2["query"].strand,
        )
        return self.contig.align((loc_ref, loc_query), align_method)

    @classmethod
    def _parse_alignment(cls, alignment: Iterable[int]):
        """
        this method follows the format of the delta file:

        >   Each of these headers is followed by a string of signed digits, one per line,
        with the final line before the next header equaling 0 (zero). Each digit
        represents the distance to the next insertion in the reference (positive int)
        or deletion in the reference (negative int), as measured in DNA bases or amino
        acids depending on the alignment data type. For example, with 'nucmer' the
        delta sequence (1, -3, 4, 0) would represent an insertion at positions 1 and 7
        in the reference sequence and an insertion at position 3 in the query sequence.
        Or with letters:

        > A = acgtagctgag$
        > B = cggtagtgag$
        > Delta = (1, -3, 4, 0)
        > A = acg.tagctgag$
        > B = .cggtag.tgag$

        >>> AlignRegion._parse_alignment([1, -3, 4])
        ([1, 8], [4])
        """
        ref: list[int] = []
        query: list[int] = []
        prev = 0
        for i in alignment:
            if i < 0:
                prev += -i
                query.append(prev)
            elif i == 0:  # pragma: no cover
                break
            else:
                prev += i
                ref.append(prev)
        return ref, query
