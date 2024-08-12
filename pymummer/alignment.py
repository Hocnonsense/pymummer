# -*- coding: utf-8 -*-
"""
 * @Date: 2024-08-11 17:37:59
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-12 11:17:10
 * @FilePath: /pymummer/pymummer/alignment.py
 * @Description:
"""
# """

import os
from typing import Iterable, Literal, overload

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .pair import Pair, align_edlib


class AlignContig2:
    def __init__(
        self,
        seqids: tuple[str, str],
        seqs: tuple[SeqRecord, SeqRecord],
    ):
        ref_id, query_id = seqids
        ref_seq, query_seq = seqs
        self.seqid2 = Pair({"ref": ref_id, "query": query_id})
        self.seq2 = Pair({"ref": ref_seq, "query": query_seq})
        self.alignregions: list[AlignRegion] = []

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
        strand: tuple[Literal[1, -1], Literal[1, -1]] = (1, 1),
        /,
        align_method: Literal["edlib"] = "edlib",
    ) -> "AlignRegion": ...
    @overload
    def align(
        self,
        locs: tuple[SimpleLocation, SimpleLocation],
        /,
        align_method: Literal["edlib"] = "edlib",
    ) -> "AlignRegion": ...
    def align(self, locs=1, /, align_method="edlib"):
        # check locs
        if isinstance(locs, int):
            loc_ref = SimpleLocation(0, len(self.seq2["ref"]), 1)
            loc_query = SimpleLocation(0, len(self.seq2["query"]), locs)
        elif isinstance(locs, tuple):
            _ref, _query = locs
            if isinstance(_ref, int):
                loc_ref = SimpleLocation(0, len(self.seq2["ref"]), _ref)
            else:
                loc_ref = _ref
            if isinstance(_query, int):
                loc_query = SimpleLocation(0, len(self.seq2["query"]), _query)
            else:
                loc_query = _query
        else:
            raise NotImplementedError(f"not implied {locs = }")  # pragma: no cover
        # check align_method and align
        if align_method == "edlib":
            ref_seq = loc_ref.extract(self.seq2["ref"].seq)
            query_seq = loc_query.extract(self.seq2["query"].seq)
            aln, mis = align_edlib(ref_seq, query_seq)
        else:
            raise NotImplementedError(
                f"not implied {align_method = }"
            )  # pragma: no cover
        # return AlignRegion
        return self.region_class((loc_ref, loc_query), aln, mis, self)

    @classmethod
    def _get_region_class(cls):
        return AlignRegion

    @property
    def region_class(self):
        return self._get_region_class()


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
        refid = self.contig.seqid2["ref"] if self.contig else None
        queryid = self.contig.seqid2["query"] if self.contig else None
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
            f" {self.seq_len}bp-{self.n_mismatch})"
        )

    @property
    def seq_len(self):
        return max(len(i) for i in self.loc2.values())

    @property
    @Pair
    def seq(self, this: "Pair.T", other: "Pair.T") -> Seq | None:
        if self.contig is not None:
            seq = self.contig.seq2[this]
            if seq is not None:
                return self.loc2[this].extract(seq.seq)
        return None  # pragma: no cover

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
        if seq is not None:
            return self._align_indel(seq, self.feat2[this].qualifiers["ins"])
        return None  # pragma: no cover

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
        if ref is None or query is None:
            return None  # pragma: no cover
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
        assert not ins, "unexpect base left "

    @property
    @Pair
    def diff2(self, this: "Pair.T", other: "Pair.T"):
        ref, query = self.seq_align[this], self.seq_align[other]
        if ref is None or query is None:
            return None  # pragma: no cover
        return self.pair2diff(ref, query)

    @property
    def diff(self):
        return self.diff2["ref"]

    @classmethod
    def mask_muts(cls, alignment: str, n_wing=5, min_mask=10):
        """
        >>> AlignRegion.mask_muts("||||||+|||-|||||", 1, 3)
        [(0, 5), (12, 16)]
        """
        last_char = ""
        start_i = 0
        align_regions: list[tuple[int, int]] = []
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
            self.loc2[i].strand == other.loc2[i].strand for i in ("ref", "query")
        ), f"Strand not the same, this: {self}, other: {other}"
        assert self.contig is not None
        loc_ref = SimpleLocation(
            min(self.loc2["ref"].start, other.loc2["ref"].start),
            max(self.loc2["ref"].end, other.loc2["ref"].end),
            self.loc2["ref"].strand,
        )
        loc_query = SimpleLocation(
            min(self.loc2["query"].start, other.loc2["query"].start),
            max(self.loc2["query"].end, other.loc2["query"].end),
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
        """
        ref: list[int] = []
        query: list[int] = []
        prev = 0
        for i in alignment:
            if i < 0:
                prev += -i
                query.append(prev)
            else:
                prev += i
                ref.append(prev)
        return ref, query
