from typing import (
    Callable,
    Generic,
    Iterable,
    Literal,
    Sequence,
    TypeVar,
    Final,
    overload,
)
from Bio.Align import PairwiseAligner

FT = TypeVar("FT")
CLS = TypeVar("CLS")


# class Pair(Generic[FT, CLS]):
class Pair(Generic[FT]):
    """
    >>> class Test:
    ...     @property
    ...     @Pair
    ...     def l(self, this: "Pair.T", other: "Pair.T"):
    ...         return other[::-1]
    ...
    >>> t = Test()
    >>> a = t.l["ref"]
    >>> assert a == "yreuq"
    >>> t.l.query
    'fer'
    >>> t.l["X"]
    Traceback (most recent call last):
        ...
    KeyError: '"X" not in switch'
    """

    T = Literal["ref", "query"]

    @overload
    def __init__(self, f: Callable[[CLS, T, T], FT]): ...

    @overload
    def __init__(self, f: dict[T, FT]): ...

    def __init__(
        self, f: Callable[[CLS, T, T], FT] | dict[T, FT], inst: CLS | None = None
    ):
        if isinstance(f, dict):
            self._d = f
        else:
            self._f = f
        self.inst = inst

    def __call__(self, inst):
        if hasattr(self, "_d"):
            raise NotImplementedError("Cannot call Pair with dict")
        self.inst = inst
        return self

    ENUM: Final[tuple[Literal["ref"], Literal["query"]]] = T.__args__  # type: ignore

    @classmethod
    def switch(cls, this: T):
        if this == cls.ENUM[0]:
            return cls.ENUM[1]
        if this == cls.ENUM[1]:
            return cls.ENUM[0]
        raise KeyError(f'"{this}" not in switch')

    def __getitem__(self, this: T):
        if self.inst is None:
            return self._d[this]
        return self._f(self.inst, this, self.switch(this))

    @property
    def ref(self):
        return self["ref"]

    @property
    def query(self):
        return self["query"]

    def __iter__(self):
        yield from self.ENUM

    def values(self):
        for i in self:
            yield self[i]

    def items(self):
        for i in self:
            yield i, self[i]


alner_extend_gap = PairwiseAligner(
    match_score=1,
    mismatch_score=-1,
    open_gap_score=-1,
    extend_gap_score=-0.1,
    mode="global",
)


def align_biopair(seq1, seq2, alner=alner_extend_gap):
    """
    compare two sequences using biopython
    return similar format that MUMmer likes:
    >>> a = "acgtagctgagacgtagctgaggtgagk"
    >>> b = "acgtagctgaggtgagck"

    >>> align_biopair("acgtagctgag", "cggtagtgag")
    ([1, -2, 5], 3)
    >>> align_biopair("acgtagctgagacgtagctgaggtgagk", "acgtagctgaggtgagck")
    ([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -17], 12)
    """
    aln = next(alner.align(seq1, seq2))
    aln1, aln2 = list(aln)  # pyright: ignore[reportArgumentType]
    return aln2delta(aln1, aln2)


def aln2delta(aln1, aln2):
    muts: list[int] = []
    ndiff = 0
    match_base = 0
    for a, b in zip(aln1, aln2):
        match_base += 1
        if a != b:
            ndiff += 1
            if a == "-":
                muts.append(-match_base)
            elif b == "-":
                muts.append(match_base)
            else:
                continue
            match_base = 0
    return muts, ndiff


def diff2delta(diff: Iterable[str]):
    """
    > A = acg.tagctgag$
    > B = .cggtag.tgag$
    >>> diff2delta(["-", "|", "|", "g+|", "|", "|", "-", "|", "|", "|", "|"])
    ([1, -3, 4], 3)
    >>> diff2delta(["|", "|+C", ])
    ([-3], 1)
    """
    muts: list[int] = []
    ndiff = 0
    match_base = 0
    for c in diff:
        if c == "-":
            muts.append(match_base + 1)
            match_base = 0
            ndiff += 1
        elif c == "|":
            match_base += 1
        elif len(c) == 1:
            match_base += 1
            ndiff += 1
        elif "+" in c:
            if c.startswith("|+"):
                muts.append(-match_base - 2)
            elif c.endswith("+|"):
                muts.append(-match_base - 1)
                match_base = 1
            else:
                raise ValueError(f"Invalid diff: {c}")
            ndiff += 1
            for i in range(len(c) - 3):
                muts.append(-1)
                ndiff += 1
    return muts, ndiff


def cigar2delta(cigar: str):
    muts: list[int] = []
    ndiff = 0
    match_base = 0
    last_str = cigar[0]
    for c in cigar[1:]:
        if c == "=":
            # equal for equal or mismatch
            match_base += int(last_str)
        elif c == "X":
            match_base += int(last_str)
            ndiff += int(last_str)
        elif c == "I":
            # positive to insert to seq2
            muts.append(match_base + 1)
            for i in range(int(last_str) - 1):
                muts.append(1)
                ndiff += 1
            match_base = 0
        elif c == "D":
            # negative to insert to seq1
            muts.append(-match_base - 1)
            for i in range(int(last_str) - 1):
                muts.append(-1)
                ndiff += 1
            match_base = 0
        else:
            # recort number
            last_str += c
            continue
        last_str = ""
    return muts, ndiff


try:
    import edlib

    def align_edlib(seq1, seq2):
        """
        compare two sequences using edlib
        return similar format that MUMmer likes:
        > A = acgtagctgag$
        > B = cggtagtgag$
        > Delta = (1, -3, 4, 0)
        > A = acg.tagctgag$
        > B = .cggtag.tgag$
              I==D===I====$
        cigar: 1I2=1D3=1I4=

        >>> align_edlib("acgtagctgag", "cggtagtgag")
        ([1, -3, 4], 3)
        """
        seq1, seq2 = str(seq1), str(seq2)
        assert edlib is not None
        ealign = edlib.align(seq1, seq2, mode="NW", task="path")
        # print(ealign)
        ndiff = int(ealign["editDistance"])
        muts = cigar2delta(ealign["cigar"])[0]
        return muts, ndiff

    IMPORT_AVAIL_EDLIB = True
except ImportError:  # pragma: no cover
    IMPORT_AVAIL_EDLIB = False


try:
    from hgvs import sequencevariant, posedit, location, edit

    DOC_DEFINE_HGVS = """Creates a genomic SequenceVariant from a VCF record and the specified alt

    https://www.hgvs.org/
    https://onlinelibrary.wiley.com/doi/10.1002/humu.22981

    Table 1. Nomenclature Definitions with Example Variant Descriptions
    | Type          | Example               | Description
    |---------------|-----------------------|----------------------------------------
    | [>][2]        | g.1318G>T             | […][1] 1 nucleotide is replaced by 1 other nucleotide
    | del(etion)	| g.3661_3706del        | … >=1 nucleotides are not present (deleted)
    | ins(ertion)	| g.7339_7340insTAGG    | … >=1 nucleotides are inserted in a sequence and where the insertion is not a copy of a sequence immediately 5'
    | [indel][3]    | g.112_117delinsTG     | … >=1 nucleotides are replaced by >=1 other nucleotides and which is not a substitution, inversion, or conversion
    | inv(ersion)	| g.495_499inv          | … more than 1 nucleotide replaces the original sequence and is the reverse-complement of the original sequence (e.g., CTCGA to TCGAG)
    | dup(location)	| g.3661_3706dup        | … a copy of >=1 nucleotides are inserted directly 3' of the original copy of that sequence
    | con(version)	| g.333_590con1844_2101 | A specific type of deletion-insertion where a range of nucleotides replacing the original sequence are a copy of a sequence from another site in the genome

    [1]: … "a change where in a specific sequence compared to the reference sequence"
    [2]: > "Substitution"
    [3]: indel "Deletion-insertion"
    """

    def hgvs_from_mut(
        start: int, ref: Sequence, alt: Sequence, seqid: str, seqtype="g"
    ):
        """
        I try to generate a hgvs record as a [vcf input](https://github.com/biocommons/hgvs/blob/83f297ddf10132f9a51b2163fa64b4204b40586f/misc/experimental/vcf-add-hgvs)
        - a vcf.model._Record has this (pip install pyvcf):
            ```python
            _Record(
                CHROM=chrom,
                POS=pos,
                REF=ref,
                ALT=alt,
                QUAL=qual,
                FILTER=filt,
                INFO=info,
                FORMA=fmt,
            )
            self.start = self.POS - 1
            #: zero-based, half-open end coordinate of ``REF``
            self.end = self.start + len(self.REF)
            ```
        - and transfer from _Record to SequenceVariant:
            ```python
            ref = r.REF
            alt = r.ALT[alt_index].sequence if r.ALT[alt_index] else ""
            start = r.POS - 1
            end = start + len(r.REF)
            ac = chrome2seqid[r.CHROM]
            ```
        """
        end = start + len(ref)
        if ref == "" and alt != "":
            # insertion
            end += 1
        else:
            start += 1
        # pfx = os.path.commonprefix([ref, alt])
        # lp = len(pfx)
        # if lp > 0:
        #    ref = ref[lp:]
        #    alt = alt[lp:]
        #    start += lp
        return sequencevariant.SequenceVariant(
            seqid,  # ac  # pyright: ignore[reportCallIssue]
            seqtype,  # type
            posedit.PosEdit(  # posedit
                location.Interval(  # pyright: ignore[reportCallIssue]
                    location.SimplePosition(start),  # pyright: ignore[reportCallIssue]
                    location.SimplePosition(end),  # pyright: ignore[reportCallIssue]
                    False,  # uncertain
                ),
                edit.NARefAlt(
                    ref or None, alt or None, False  # pyright: ignore[reportCallIssue]
                ),
            ),
        )

    IMPORT_AVAIL_HGVS = True
except ImportError:  # pragma: no cover
    IMPORT_AVAIL_HGVS = False
