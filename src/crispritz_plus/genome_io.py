"""
FASTA genome reading and writing for CRISPRitz-plus.
 
Provides :class:`GenomeReader`, which loads a single-record FASTA file into
memory and supports in-place SNP insertion and indel-window extraction for
the variant-enrichment pipeline, and :class:`GenomeWriter`, which writes a
header / sequence pair back out as FASTA.
 
Two parallel sequence buffers are maintained by :class:`GenomeReader`:
 
``sequence``
    The pristine, upper-cased reference sequence as read from disk.
``sequence_enr``
    A mutable copy ("enriched") into which SNP alleles are written via
    :meth:`GenomeReader.insert_snp` without disturbing the reference.
 
Module-level constants
----------------------
INDELOFFSET : int
    Number of reference bases retained upstream **and** downstream of an
    indel position when extracting its local window (see
    :meth:`GenomeReader.insert_indel`).
"""

from typing import List

import os


from .crispritz_errors import GenomeReaderError, GenomeWriterError
from .exception_handlers import exception_handler
from .enrichment.variants import IndelPair


#: Number of reference bases retained both upstream and downstream of an
#: indel position when extracting its local window.
# define indel upstream and downstream offset
INDELOFFSET = 50


# ==============================================================================
# Public classes
# ==============================================================================


class GenomeReader:
    """Load a single-record FASTA and support variant enrichment in memory.
 
    On :meth:`read`, the file is parsed into a header string and a list of
    single-character bases.  Two buffers are kept: an immutable reference
    (:attr:`sequence`) and a mutable enriched copy (:attr:`sequence_enr`)
    into which SNPs are written.
 
    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file to read.
    debug : bool
        When *True*, errors propagate with a full traceback instead of a
        formatted user-facing message.
 
    Attributes
    ----------
    header : str
        FASTA header line (without the leading ``'>'``), populated by
        :meth:`read`.
    sequence : List[str]
        The pristine reference sequence as a list of bases.
    sequence_enr : List[str]
        The mutable, SNP-enriched copy of the reference sequence.
    fname : str
        The source FASTA path.
    """

    def __init__(self, fasta_path: str, debug: bool):
        self._debug = debug  # store debug flag
        self._fasta_path = fasta_path
        self._sequence: List[str] = []
        self._sequence_enr: List[str] = []
        self._header: str = ""

    def __repr__(self) -> str:
        seq_preview = f"[{len(self._sequence)} bases]" if self._sequence else "None"
        return (
            f"GenomeReader(fasta_path={self._fasta_path!r}, "
            f"header={self._header!r}, "
            f"sequence={seq_preview}, "
            f"debug={self._debug!r})"
        )

    def __str__(self) -> str:
        if not self._header or not self._sequence:
            status = "not read"
            details = ""
        else:
            seq_len = len(self._sequence)
            seq_preview = "".join(self._sequence[:50])
            if seq_len > 50:
                seq_preview += f"... ({seq_len - 50} more bases)"
            status = "read"
            details = f"\n  Header: {self._header}\n  Length: {seq_len} bases\n  Preview: {seq_preview}"
        return f"GenomeReader: {self._fasta_path} ({status}){details}"

    def _extract_header(self, header: str) -> None:
        """Store the parsed header, validating that it is non-empty.
 
        Parameters
        ----------
        header : str
            The header text (already stripped of the leading ``'>'`` and
            surrounding whitespace).
 
        Raises
        ------
        GenomeReaderError
            If *header* is empty.
        """
        self._header = header
        if not self._header:
            exception_handler(
                GenomeReaderError,
                f"FASTA header is empty: {self._fasta_path}",
                os.EX_IOERR,
                self._debug,
            )

    def _extract_sequence(self, lines: List[str]) -> None:
        """Parse, upper-case, and store the sequence body.
 
        Filters out any residual header lines, concatenates the remaining
        lines, upper-cases the result, and stores it as a list of bases in
        both :attr:`sequence` and the enriched copy :attr:`sequence_enr`.
 
        Parameters
        ----------
        lines : List[str]
            The non-header lines of the FASTA file.
 
        Raises
        ------
        GenomeReaderError
            If no sequence data is present.
        """
        sequence_lines = [line for line in lines if not line.startswith(">")]
        if not sequence_lines:
            exception_handler(
                GenomeReaderError,
                f"FASTA file contains no sequence data: {self._fasta_path}",
                os.EX_IOERR,
                self._debug,
            )
        # join and convert to list of characters
        self._sequence = list("".join(sequence_lines).upper())
        self._sequence_enr = list(self._sequence)

    def read(self) -> None:
        """Read and parse the FASTA file into header and sequence buffers.
 
        Strips blank lines, requires the first line to be a header
        (``'>'``-prefixed), and populates :attr:`header`, :attr:`sequence`,
        and :attr:`sequence_enr`.
 
        Returns
        -------
        None
            Populates the reader's internal state as a side-effect.
 
        Raises
        ------
        GenomeReaderError
            If the file is empty, does not start with a header line, has an
            empty header, contains no sequence data, or cannot be read.
        """
        try:
            with open(self._fasta_path, mode="r") as fin:
                lines = [line.rstrip("\n\r") for line in fin if line.strip()]
                if not lines:
                    exception_handler(
                        GenomeReaderError,
                        f"FASTA file is empty: {self._fasta_path}",
                        os.EX_IOERR,
                        self._debug,
                    )
                if not lines[0].startswith(">"):
                    exception_handler(
                        GenomeReaderError,
                        f"FASTA file must start with a header line: {self._fasta_path}",
                        os.EX_IOERR,
                        self._debug,
                    )
                # extract header (remove '>' and strip whitespace)
                self._extract_header(lines[0].lstrip(">").strip())
                # extract sequence (all non-header lines)
                self._extract_sequence(lines[1:])
        except (IOError, Exception) as e:
            exception_handler(
                GenomeReaderError,
                f"Failed reading FASTA: {self._fasta_path}",
                os.EX_IOERR,
                self._debug,
                e,
            )

    def insert_snp(self, iupac_nt: str, pos: int) -> None:
        """Write a SNP allele into the enriched sequence at *pos*.
 
        Overwrites position *pos* of :attr:`sequence_enr` with *iupac_nt*,
        leaving the pristine :attr:`sequence` untouched.
 
        Parameters
        ----------
        iupac_nt : str
            The IUPAC nucleotide symbol to write (a single character).
        pos : int
            Zero-based position in the enriched sequence to overwrite.
 
        Returns
        -------
        None
 
        Raises
        ------
        GenomeReaderError
            If the sequence has not yet been initialised (i.e. :meth:`read`
            has not been called).
        """
        if not self._sequence_enr:
            exception_handler(
                GenomeReaderError,
                "Sequence not initialized, impossible to insert SNPs",
                os.EX_DATAERR,
                self._debug,
            )
        self._sequence_enr[pos] = iupac_nt  # add snp as iupac to sequence

    def insert_indel(self, indel: str, pos: int, offset: int) -> IndelPair:
        """Extract the reference window around *pos* and build its indel variant.
 
        Slices a window of :attr:`sequence` spanning :data:`INDELOFFSET` bases
        upstream of *pos* through ``INDELOFFSET + offset`` bases downstream,
        then constructs the corresponding indel sequence by replacing the
        ``offset``-length reference segment (immediately after the upstream
        flank) with the supplied *indel* bases.
 
        Parameters
        ----------
        indel : str
            The indel sequence to splice in, replacing the reference segment.
        pos : int
            Zero-based reference position around which to centre the window.
        offset : int
            Length of the reference segment that the indel replaces, used to
            size the downstream extent of the extracted window.
 
        Returns
        -------
        IndelPair
            A pair holding the extracted reference window (``refseq``) and the
            corresponding indel-applied window (``indelseq``).
 
        Notes
        -----
        Both flanks are taken from the pristine :attr:`sequence`, not the
        enriched copy.  The reference and indel windows share the same
        upstream and downstream flanks; only the central segment differs.
        """
        refseq = self._sequence[pos - INDELOFFSET : pos + INDELOFFSET + offset]
        indelseq = refseq[:INDELOFFSET] + list(indel) + refseq[INDELOFFSET + offset :]
        return IndelPair(refseq=refseq, indelseq=indelseq)

    def to_string(self) -> str:
        """Return the enriched sequence as a single string.
 
        Returns
        -------
        str
            The concatenation of :attr:`sequence_enr`.
 
        Raises
        ------
        GenomeReaderError
            If the sequence has not yet been initialised.
        """
        if not self._sequence_enr:
            exception_handler(
                GenomeReaderError,
                f"Sequence not initialized, impossible to retrieve sequence as {str.__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        return "".join(self._sequence_enr)

    @property
    def header(self) -> str:
        """str: The FASTA header (without the leading ``'>'``).
 
        Raises
        ------
        AssertionError
            If accessed before :meth:`read` populates the header.
        """
        assert self._header  # if used, always initialized
        return self._header

    @property
    def sequence(self) -> List[str]:
        """List[str]: The pristine reference sequence as a list of bases.
 
        Raises
        ------
        AssertionError
            If accessed before :meth:`read` populates the sequence.
        """
        assert self._sequence  # if used, always initialized
        return self._sequence

    @property
    def sequence_enr(self) -> List[str]:
        """List[str]: The mutable, SNP-enriched copy of the sequence.
 
        Raises
        ------
        AssertionError
            If accessed before :meth:`read` populates the sequence.
        """
        assert self._sequence_enr  # if used, always initialized
        return self._sequence_enr

    @property
    def fname(self) -> str:
        """str: The source FASTA path."""
        return self._fasta_path


class GenomeWriter:
    """Write a header / sequence pair to disk as a single-record FASTA.
 
    Parameters
    ----------
    outfile : str
        Destination path for the FASTA output.
    debug : bool
        When *True*, write errors propagate with a full traceback instead of
        a formatted user-facing message.
    """

    def __init__(self, outfile: str, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._outfile = outfile

    def __repr__(self) -> str:
        return f"GenomeWriter(outfile={self._outfile!r}, debug={self._debug!r})"

    def __str__(self) -> str:
        return f"GenomeWriter: {self._outfile}"

    def write(self, header: str, sequence_list: List[str]) -> None:
        """Write *header* and *sequence_list* to *outfile* in FASTA format.
 
        Emits a ``'>'``-prefixed header line followed by the joined sequence
        on a single line, each terminated by a newline.
 
        Parameters
        ----------
        header : str
            Header text to write (the ``'>'`` prefix is added automatically).
        sequence_list : List[str]
            The sequence as a list of single-character strings; joined into
            one line before writing.
 
        Returns
        -------
        None
 
        Raises
        ------
        GenomeWriterError
            If the file cannot be written.
        """
        sequence = "".join(sequence_list)
        try:
            with open(self._outfile, mode="w") as fout:
                fout.write(f">{header}\n")  # write header
                fout.write(f"{sequence}\n")  # write sequence with newline
        except (IOError, Exception) as e:
            exception_handler(
                GenomeWriterError,
                f"Failed writing FASTA: {self._outfile}",
                os.EX_IOERR,
                self._debug,
                e,
            )
