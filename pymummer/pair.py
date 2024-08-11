from typing import Callable, Generic, Literal, TypeVar, Final, overload


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
    ENUM: Final[tuple[Literal["ref"], Literal["query"]]] = T.__args__  # type: ignore

    @classmethod
    def __iter__(cls):
        yield from cls.ENUM

    @overload
    def __init__(self, f: Callable[[CLS, T, T], FT]): ...

    @overload
    def __init__(self, f: dict[T, FT]): ...

    def __init__(
        self, f: Callable[[CLS, T, T], FT] | dict[T, FT], cls: CLS | None = None
    ):
        if isinstance(f, dict):
            self._d = f
        else:
            self._f = f
        self.cls = cls

    def __call__(self, cls: CLS):
        self.cls = cls
        return self

    @classmethod
    def switch(cls, this: T):
        if this == cls.ENUM[0]:
            return cls.ENUM[1]
        if this == cls.ENUM[1]:
            return cls.ENUM[0]
        raise KeyError(f'"{this}" not in switch')

    def __getitem__(self, this: T):
        if self.cls is None:
            return self._d[this]
        return self._f(self.cls, this, self.switch(this))

    @property
    def ref(self):
        return self["ref"]

    @property
    def query(self):
        return self["query"]


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

        >>> analyse_edlib("acgtagctgag", "cggtagtgag")
        ([1, -3, 4], 3)
        """
        seq1, seq2 = str(seq1), str(seq2)
        ealign = edlib.align(seq1, seq2, mode="NW", task="path")
        print(ealign)
        ndiff = int(ealign["editDistance"])
        cigar: str = ealign["cigar"]
        muts: list[int] = []
        match_base = 0
        last_str = cigar[0]
        for c in cigar[1:]:
            if c == "=":
                # equal for equal or mismatch
                match_base += int(last_str)
            elif c == "X":
                match_base += int(last_str)
            elif c == "I":
                # positive to insert to seq2
                muts.append(match_base + 1)
                for i in range(int(last_str) - 1):
                    muts.append(1)
                match_base = 0
            elif c == "D":
                # negative to insert to seq1
                muts.append(-match_base - 1)
                for i in range(int(last_str) - 1):
                    muts.append(-1)
                match_base = 0
            else:
                # recort number
                last_str += c
                continue
            last_str = ""
        return muts, ndiff

except ImportError:
    edlib = None
