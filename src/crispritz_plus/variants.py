"""
Variant carrier types for the CRISPRitz-plus enrichment pipeline.

Provides lightweight, mostly-immutable representations of the SNPs and
indels parsed from a VCF record, their per-record collections, and the
supporting value objects used when materialising indels into the synthetic
("fake") FASTA used during enrichment.

Class overview
--------------
Snp / Snps
    A single SNP allele and a per-record collection of SNP alleles.
Indel / Indels
    A single indel allele (with a symbolic-allele check) and a per-record
    collection of indel alleles (with symbolic filtering).
IndelInfo
    Coordinate record locating one indel within the synthetic FASTA.
IndelPair
    A reference / indel sequence pair produced when an indel is inserted.
IndelsSet
    Accumulator that lays out indel sequences end-to-end in the synthetic
    FASTA and tracks their coordinates.

Coordinate convention
---------------------
All variant positions are **0-based**, matching the convention used by the
rest of the enrichment code.
"""

from dataclasses import dataclass
from typing import List, Iterator, Optional


# ==============================================================================
# Data classes
# ==============================================================================


@dataclass
class IndelInfo:
    """Coordinate record locating one indel within the synthetic FASTA.

    Attributes
    ----------
    idx : int
        Identifier assigned to the indel by :class:`IndelsSet`.
    start : int
        Start offset of the indel sequence within the synthetic FASTA.
    stop : int
        Stop offset (``start + len(indel_seq)``) within the synthetic FASTA.
    """

    idx: int
    start: int
    stop: int


@dataclass
class IndelPair:
    """Reference / indel sequence pair produced when an indel is inserted.

    Attributes
    ----------
    refseq : List[str]
        The reference window (as a list of single-character strings)
        surrounding the indel position.
    indelseq : List[str]
        The same window with the indel applied.
    """

    refseq: List[str]
    indelseq: List[str]


# ==============================================================================
# Public Classes
# ==============================================================================


class Snp:
    """Represent a single SNP allele from a VCF record. This stores its genomic
    position, reference and alternate bases, and the genotype index pointing to
    the allele.

    The object is immutable after creation and is used as a lightweight carrier
    when building SNP collections for enrichment and annotation.
    """

    def __init__(self, pos: int, ref: str, alt: str, gtidx: int) -> None:
        # NOTE: pos are 0-based
        self._pos = pos  # snp position
        self._ref = ref  # ref allele
        self._alt = alt  # alt allele
        self._gtidx = gtidx  # genotype index for variant

    def __repr__(self) -> str:
        return (
            f"<Snp object; pos={self._pos}, ref={self._ref!r}, "
            f"alt={self._alt!r}, gtidx={self._gtidx}>"
        )

    @property
    def pos(self) -> int:
        """int: Zero-based genomic position of the SNP."""
        return self._pos

    @property
    def ref(self) -> str:
        """str: Reference base from the VCF record."""
        return self._ref

    @property
    def alt(self) -> str:
        """str: Alternate base for this allele."""
        return self._alt

    @property
    def gtidx(self) -> int:
        """int: Genotype index referring to this allele."""
        return self._gtidx


class Snps:
    """Collect multiple SNP alleles that originate from the same VCF record.
    This groups individual `Snp` objects and exposes convenience accessors over
    the set.

    The container behaves like a lightweight sequence with iteration, truthiness,
    and helper methods to retrieve shared position, reference base and allele
    metadata across all SNPs.
    """

    def __init__(self, items: Optional[List[Snp]] = None) -> None:
        self._items: List[Snp] = items if items is not None else []

    def __repr__(self) -> str:
        return f"Snps(n={len(self._items)}, items={self._items!r})"

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[Snp]:
        return iter(self._items)

    def __bool__(self) -> bool:
        return bool(self._items)

    def add(self, snp: Snp) -> None:
        """Append a SNP allele to the collection.

        Parameters
        ----------
        snp : Snp
            The allele to add.

        Raises
        ------
        TypeError
            If *snp* is not a :class:`Snp` instance.
        """
        if not isinstance(snp, Snp):
            raise TypeError(f"Snps.add expects a Snp, got {type(snp).__name__}")
        self._items.append(snp)

    @property
    def items(self) -> List[Snp]:
        """List[Snp]: A shallow copy of the stored alleles.

        A copy is returned so that external code cannot mutate the internal
        list in place.
        """
        # expose a copy to avoid accidental external mutation
        return list(self._items)

    @property
    def alts(self) -> List[str]:
        """List[str]: The alternate base of every stored allele, in order."""
        return [s.alt for s in self._items]

    @property
    def gtidxs(self) -> List[int]:
        """List[int]: The genotype index of every stored allele, in order."""
        return [s.gtidx for s in self._items]

    @property
    def pos(self) -> int:
        """int: The shared 0-based position of the alleles in this record.

        Raises
        ------
        AssertionError
            If the collection is empty.
        """
        # all SNPs from one record share the same position
        assert self._items
        return self._items[0].pos

    @property
    def ref(self) -> str:
        """str: The shared reference base of the alleles in this record.

        Raises
        ------
        AssertionError
            If the collection is empty.
        """
        assert self._items
        return self._items[0].ref


class Indel:
    """Represent a single indel allele from a VCF record. This stores its genomic
    position, reference and alternate sequences, and the genotype index for the
    allele.

    The object also exposes a convenience check for symbolic alleles so that
    callers can consistently skip non-literal representations.
    """

    def __init__(self, pos: int, ref: str, alt: str, gtidx: int) -> None:
        # NOTE: pos are 0-based
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._gtidx = gtidx

    def __repr__(self) -> str:
        return f"<Indel object; pos={self._pos}, ref={self._ref!r}, alt={self._alt!r}, gtidx={self._gtidx}>"

    @property
    def pos(self) -> int:
        """int: Zero-based genomic position of the indel."""
        return self._pos

    @property
    def ref(self) -> str:
        """str: Reference allele sequence from the VCF record."""
        return self._ref

    @property
    def alt(self) -> str:
        """str: Alternate allele sequence from the VCF record."""
        return self._alt

    @property
    def gtidx(self) -> int:
        """int: Genotype index referring to this allele."""
        return self._gtidx

    def is_symbolic(self) -> bool:
        """Return whether either allele is a symbolic (``<...>``) allele.

        Symbolic alleles (e.g. ``<DEL>``, ``<DUP>``) do not carry a literal
        sequence and are skipped by downstream enrichment.  Centralising the
        check here keeps that rule in one place.

        Returns
        -------
        bool
            ``True`` if ``'<'`` appears in either the reference or alternate
            allele, ``False`` otherwise.
        """
        # skipped <...> alleles; this helps keep that rule centralized
        return "<" in self._ref or "<" in self._alt


class Indels:
    """Collect multiple indel alleles that originate from the same VCF record.
    This groups individual `Indel` objects and provides helpers to access shared
    coordinates and filter alleles.

    The container behaves like a lightweight sequence with iteration, truthiness
    and convenience methods for retrieving positions, reference bases and
    non-symbolic subsets of the stored indels.
    """

    def __init__(self, items: Optional[List[Indel]] = None) -> None:
        self._items: List[Indel] = items if items is not None else []

    def __repr__(self) -> str:
        return f"Indels(n={len(self._items)}, items={self._items!r})"

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[Indel]:
        return iter(self._items)

    def __bool__(self) -> bool:
        return bool(self._items)

    def add(self, indel: Indel) -> None:
        """Append an indel allele to the collection.

        Parameters
        ----------
        indel : Indel
            The allele to add.

        Raises
        ------
        TypeError
            If *indel* is not an :class:`Indel` instance.
        """
        if not isinstance(indel, Indel):
            raise TypeError(f"Indels.add expects an Indel, got {type(indel).__name__}")
        self._items.append(indel)

    @property
    def items(self) -> List[Indel]:
        """List[Indel]: A shallow copy of the stored alleles."""
        return list(self._items)

    def alts(self) -> List[str]:
        """Return the alternate sequence of every stored allele, in order.

        Returns
        -------
        List[str]
            One alternate-allele sequence per stored indel.
        """
        return [i.alt for i in self._items]

    def gtidxs(self) -> List[int]:
        """Return the genotype index of every stored allele, in order.

        Returns
        -------
        List[int]
            One genotype index per stored indel.
        """
        return [i.gtidx for i in self._items]

    def pos(self) -> int:
        """Return the shared 0-based position of the alleles in this record.

        Returns
        -------
        int
            The position of the first stored indel (all share the same one).

        Raises
        ------
        AssertionError
            If the collection is empty.
        """
        assert self._items
        return self._items[0].pos

    def ref(self) -> str:
        """Return the shared reference sequence of the alleles in this record.

        Returns
        -------
        str
            The reference allele of the first stored indel.

        Raises
        ------
        AssertionError
            If the collection is empty.
        """
        assert self._items
        return self._items[0].ref

    def non_symbolic(self) -> "Indels":
        """Return a new collection with all symbolic alleles removed.

        Filters out any allele for which :meth:`Indel.is_symbolic` is
        ``True``, leaving only literal-sequence indels.

        Returns
        -------
        Indels
            A new :class:`Indels` containing only the non-symbolic alleles.
        """
        return Indels([i for i in self._items if not i.is_symbolic()])


class IndelsSet:
    """Accumulator that lays out indel sequences end-to-end in a synthetic FASTA.

    Each pushed indel sequence is appended to an internal list and assigned a
    ``[start, stop]`` coordinate span within the synthetic ("fake") FASTA, with
    a one-base gap left between consecutive indels (``start`` of the next entry
    is ``stop + 1``).

    Parameters
    ----------
    debug : bool
        Stored debug flag for downstream error handling.
    """

    def __init__(self, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._start = 0  # indel start index in fake fasta
        self._indel_seqs: List[List[str]] = []  # indel sequences
        self._i = 1  # indel index

    def push(self, indel_seq: List[str]) -> IndelInfo:
        """Append an indel sequence and return its synthetic-FASTA coordinates.

        Stores *indel_seq*, computes its ``[start, stop]`` span at the current
        write cursor, advances the cursor past the sequence (leaving a
        one-base gap), increments the indel counter, and returns the
        coordinate record.

        Parameters
        ----------
        indel_seq : List[str]
            The indel sequence as a list of single-character strings.

        Returns
        -------
        IndelInfo
            The coordinate record for the pushed sequence.

        Notes
        -----
        The internal counter ``self._i`` is incremented **before** the
        :class:`IndelInfo` is constructed, so the returned ``idx`` reflects
        the post-increment value.
        """
        self._indel_seqs.append(indel_seq)  # push indel sequence in list
        start_i, stop_i = self._start, self._start + len(
            indel_seq
        )  # compute indel stop position in fake fasta
        self._start = stop_i + 1  # update start position
        self._i += 1  # update indel id
        return IndelInfo(idx=self._i, start=start_i, stop=stop_i)

    @property
    def start_i(self) -> int:
        """int: The current write cursor (start offset for the next push)."""
        return self._start

    @property
    def sequences(self) -> List[List[str]]:
        """List[List[str]]: All indel sequences pushed so far, in order."""
        return self._indel_seqs
