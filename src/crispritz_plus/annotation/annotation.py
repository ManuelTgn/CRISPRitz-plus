"""
Annotation module for CRISPRitz-plus.

Provides utilities to annotate CRISPR off-target search results (produced by
the ``search`` command) with genomic feature tracks supplied as BED files.
Each annotation track is sorted, block-gzip-compressed, and Tabix-indexed on
first use so that overlapping genomic features can be resolved in O(log N)
time per query via ``pysam.TabixFile``.

Pipeline overview
-----------------
1. **Track preparation** - for each BED file: sort by ``(chrom, start)`` if
   needed, bgzip-compress, and write a Tabix index (``.tbi``).  Already
   indexed tracks are passed through unchanged.  Multiple tracks are prepared
   in parallel via :class:`~concurrent.futures.ThreadPoolExecutor`.
2. **Stream annotation** - read the targets TSV line-by-line, derive each
   off-target site's 0-based half-open genomic interval, query every indexed
   track, and write one augmented TSV row per site.

Public API
----------
annotate_results
    Top-level entry point for the ``annotate-results`` CLI subcommand.

Module-level constants
----------------------
SEARCH_OUTPUT_HEADER : List[str]
    Canonical column order produced by the ``search`` command and expected
    as the first :data:`_CORE_COLUMNS` columns of every input file.
_CORE_COLUMNS : int
    Number of mandatory leading columns validated against
    :data:`SEARCH_OUTPUT_HEADER`.
_NO_OVERLAP : str
    Sentinel written to an annotation cell when no BED feature overlaps the
    off-target site (``"NA"``).
_OUTPUT_SUFFIX : str
    Extension appended to the stem of the input file to form the output
    path (``".annotated.tsv"``).
_BED_SKIP_PREFIXES : Tuple[str, ...]
    Line-prefix strings that identify non-data lines in a BED file
    (comments, ``track`` directives, and ``browser`` directives).
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

import os
import pysam


from ..exception_handlers import exception_handler
from ..utils import rename_files, remove_file
from ..verbosity import VERBOSITY_LVL, print_verbosity

from .crispritz_annotation_error import CrispritzAnnotationError


# ==============================================================================
# Module-level constants
# ==============================================================================

#: Canonical column names output by the ``search`` command, in order.
#: The first :data:`_CORE_COLUMNS` entries are validated against every input
#: file header before annotation begins.
SEARCH_OUTPUT_HEADER: List[str] = [
    "chrom",
    "pos",
    "strand",
    "grna",
    "spacer",
    "mismatches",
    "bulge_type",
    "bulge_dna",
    "bulge_rna",
    "cfd_score",
]

#: Number of leading columns in :data:`SEARCH_OUTPUT_HEADER` that **must**
#: be present and in order for a file to be accepted as a valid targets TSV.
_CORE_COLUMNS: int = 9

#: Value written to an annotation cell when no BED feature overlaps the
#: off-target site.
_NO_OVERLAP: str = "NA"

#: Suffix appended to the targets-file stem to form the output filename.
_OUTPUT_SUFFIX: str = ".annotated.tsv"

#: BED line prefixes that denote non-data lines (comment, track directive,
#: browser directive); these are skipped silently during parsing.
_BED_SKIP_PREFIXES: Tuple[str, ...] = ("#", "track", "browser")


# ==============================================================================
# Data classes
# ==============================================================================


@dataclass
class BedRecord:
    """Lightweight in-memory representation of a single BED record.

    Used exclusively during the sort-then-compress preparation phase.  The
    full original line is preserved verbatim so it can be written back after
    sorting without any field reconstruction.

    Attributes
    ----------
    chrom : str
        Chromosome / sequence name (BED column 1).
    pos : int
        0-based start position (BED column 2), used as the numeric sort key.
    line : str
        Raw, tab-separated BED line stripped of its trailing newline.
    """

    chrom: str
    pos: int
    line: str


@dataclass
class OffTarget:
    """Genomic coordinates of a single off-target site.

    Holds the half-open 0-based interval ``[start, end)`` as expected by
    ``pysam.TabixFile.fetch``.

    Attributes
    ----------
    chrom : str
        Chromosome / sequence name.
    start : int
        0-based inclusive start position.
    end : int
        0-based exclusive end position.
    strand : str
        Strand of the off-target site (``'+'`` for forward, ``'-'`` for
        reverse).
    """

    chrom: str
    start: int
    end: int
    strand: str


# ==============================================================================
# Internal helpers - output path
# ==============================================================================


def _create_targets_ann(targets_file: str, outdir: str) -> str:
    """Build the absolute output path for the annotated targets file.

    The output filename is derived from *targets_file* by stripping the last
    filename extension and appending :data:`_OUTPUT_SUFFIX`
    (``'.annotated.tsv'``).

    Parameters
    ----------
    targets_file : str
        Path to the input search-results TSV.
    outdir : str
        Directory in which the annotated file will be written.  Resolved to
        an absolute path before joining.

    Returns
    -------
    str
        Absolute path of the intended output file.
    """
    stem = os.path.basename(targets_file).rsplit(".", 1)[0]
    return os.path.join(os.path.abspath(outdir), f"{stem}{_OUTPUT_SUFFIX}")


# ==============================================================================
# Internal helpers - header validation
# ==============================================================================


def _validate_targets_header(
    header_fields: Sequence[str], debug: bool
) -> Dict[str, int]:
    """Validate the targets TSV header and return a column-index mapping.

    Checks that the first :data:`_CORE_COLUMNS` columns match
    :data:`SEARCH_OUTPUT_HEADER` exactly.  Raises
    :class:`~.crispritz_annotation_error.CrispritzAnnotationError` via
    :func:`~crispritz_plus.exception_handlers.exception_handler` when the
    schema is unexpected.

    Parameters
    ----------
    header_fields : Sequence[str]
        Column names parsed from the TSV header line.
    debug : bool
        When *True*, the exception propagates with a full traceback instead
        of being converted to a user-facing message.

    Returns
    -------
    Dict[str, int]
        Mapping of ``column_name -> 0-based column index`` for every column
        present in the header (both core and any extra columns).

    Raises
    ------
    CrispritzAnnotationError
        If the first :data:`_CORE_COLUMNS` columns do not match the
        expected schema.
    """
    expected_core = SEARCH_OUTPUT_HEADER[:_CORE_COLUMNS]
    if list(header_fields[:_CORE_COLUMNS]) != expected_core:
        exception_handler(
            CrispritzAnnotationError,
            (
                "Input file does not look like a 'search' targets table. "
                f"Expected the first {_CORE_COLUMNS} columns to be "
                f"{expected_core}, but found "
                f"{list(header_fields[:_CORE_COLUMNS])}. The annotate-results "
                "command only accepts the TSV produced by the search command"
            ),
            os.EX_DATAERR,
            debug,
        )
    return {name: idx for idx, name in enumerate(header_fields)}


def _resolve_annotation_names(
    annotation_files: List[str],
    annotation_names: Optional[List[str]],
    debug: bool,
) -> List[str]:
    """Resolve the output column headers for the annotation tracks.

    Parameters
    ----------
    annotation_files : List[str]
        The BED tracks supplied on the command line.
    annotation_names : Optional[List[str]]
        Optional explicit names, one per track.  When ``None`` the columns
        default to ``annotation1`` … ``annotationN``.
    debug : bool
        When *True*, a length-mismatch error propagates with a full
        traceback.

    Returns
    -------
    List[str]
        One column name per annotation track, in the same order as
        *annotation_files*.

    Raises
    ------
    CrispritzAnnotationError
        If *annotation_names* is provided but its length differs from that
        of *annotation_files*.
    """
    if annotation_names is None:  # fallback: annotation1 ... annotationN
        return [f"annotation{i}" for i in range(1, len(annotation_files) + 1)]
    if len(annotation_names) != len(annotation_files):
        exception_handler(
            CrispritzAnnotationError,
            (
                f"Number of annotation names ({len(annotation_names)}) does not "
                f"match the number of annotation files ({len(annotation_files)})"
            ),
            os.EX_USAGE,
            debug,
        )
    return annotation_names


# ==============================================================================
# Internal helpers - BED track preparation
# ==============================================================================


def _is_bgzipped(path: str) -> bool:
    """Return ``True`` if *path* appears to be bgzip/gzip-compressed.

    Detection is purely extension-based (``'.gz'`` suffix); no magic-byte
    inspection is performed.

    Parameters
    ----------
    path : str
        File path to inspect.

    Returns
    -------
    bool
        ``True`` when *path* ends with ``'.gz'``, ``False`` otherwise.
    """
    return path.endswith(".gz")


def _read_bed(bed_path: str, debug: bool) -> List[BedRecord]:
    """Read a plain-text BED file into a list of :class:`BedRecord` objects.

    Lines starting with ``'#'``, ``'track'``, or ``'browser'`` and blank
    lines are skipped silently.  Records with fewer than four tab-separated
    fields are also skipped (BED column 4 carries the feature label used
    during annotation and is therefore required).

    Parameters
    ----------
    bed_path : str
        Path to the plain-text BED file to read.
    debug : bool
        When *True*, I/O or parse errors propagate with a full traceback.

    Returns
    -------
    List[BedRecord]
        Parsed records in file order.

    Raises
    ------
    CrispritzAnnotationError
        On any I/O or parsing error.
    """
    try:
        records: List[BedRecord] = []
        with open(bed_path, mode="rt") as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith(_BED_SKIP_PREFIXES):
                    continue
                fields = line.split("\t")
                if len(fields) < 4:  # assume annotation on 4th column
                    continue  # not a valid BED interval record
                records.append(
                    BedRecord(chrom=fields[0], pos=int(fields[1]), line=line)
                )
    except (IOError, Exception) as e:
        exception_handler(
            CrispritzAnnotationError,
            f"An error occurred while reading: {bed_path}",
            os.EX_IOERR,
            debug,
            e,
        )
    return records


def _sort_bed(bed_path: str, debug: bool) -> str:
    """Sort a BED file by chromosome then by start position.

    Reads *bed_path* into memory via :func:`_read_bed`, sorts
    lexicographically on ``(chrom, pos)``, and writes the result to a
    sibling file named ``{bed_path}.sorted.bed``.

    Parameters
    ----------
    bed_path : str
        Path to the unsorted, plain-text BED file.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    str
        Path to the sorted plain-text BED file (``{bed_path}.sorted.bed``).

    Raises
    ------
    CrispritzAnnotationError
        On I/O or sorting failure.
    """
    try:
        records = _read_bed(bed_path, debug)
        records.sort(key=lambda rec: (rec.chrom, rec.pos))
        sorted_plain = f"{bed_path}.sorted.bed"
        with open(sorted_plain, "w") as fout:
            fout.writelines(f"{rec.line}\n" for rec in records)
    except (IOError, Exception) as e:
        exception_handler(
            CrispritzAnnotationError,
            f"Failed sorting BED {bed_path}",
            os.EX_IOERR,
            debug,
            e,
        )
    return sorted_plain


def _compress_and_index(sorted_plain: str, bed_path: str) -> str:
    """Bgzip-compress and Tabix-index a sorted plain-text BED file.

    Steps performed:

    1. Compress *sorted_plain* to ``{sorted_plain}.gz`` via
       ``pysam.tabix_compress``.
    2. Rename the compressed file to ``{bed_path}.gz``.
    3. Write a Tabix index (``{bed_path}.gz.tbi``) via
       ``pysam.tabix_index`` using the ``bed`` preset.
    4. Delete the intermediate plain-text file *sorted_plain*.

    Parameters
    ----------
    sorted_plain : str
        Path to the sorted, uncompressed BED file (typically produced by
        :func:`_sort_bed`).
    bed_path : str
        Original BED path used to derive the final ``.gz`` destination.

    Returns
    -------
    str
        Path to the bgzip-compressed, Tabix-indexed file
        (``{bed_path}.gz``).

    Notes
    -----
    The ``.tbi`` index is written alongside the ``.gz`` file by
    ``pysam.tabix_index`` following the conventional
    ``{gz_path}.tbi`` naming scheme and is not explicitly returned.
    """
    gz_path, gz_path_ = f"{sorted_plain}.gz", f"{bed_path}.gz"
    pysam.tabix_compress(sorted_plain, gz_path, force=True)
    rename_files(gz_path, gz_path_)
    pysam.tabix_index(gz_path_, preset="bed", force=True)  # writes gz_path.tbi
    remove_file(sorted_plain)
    return gz_path_


def _prepare_bed_track(bed_path: str, verbosity: int, debug: bool) -> str:
    """Ensure a BED track is bgzip-compressed and Tabix-indexed.

    Returns *bed_path* unchanged when it is already bgzip-compressed
    (``*.gz``) **and** the corresponding ``.tbi`` index exists on disk.
    Otherwise, sorts (if the file is not already compressed), compresses,
    and indexes it in-place.

    Parameters
    ----------
    bed_path : str
        Path to a BED annotation track.  May be plain-text or
        bgzip-compressed.
    verbosity : int
        Controls progress message output via
        :func:`~crispritz_plus.verbosity.print_verbosity`.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    str
        Path to the bgzip-compressed, Tabix-ready BED file.

    Raises
    ------
    CrispritzAnnotationError
        If the track cannot be prepared due to a malformed BED record, an
        I/O error, or an indexing failure.
    """
    if _is_bgzipped(bed_path) and os.path.isfile(f"{bed_path}.tbi"):
        print_verbosity(
            f"Track already indexed, skipping: {bed_path}",
            verbosity,
            VERBOSITY_LVL[2],
        )
        return bed_path  # already random-access ready; nothing to do
    print_verbosity(
        f"Preparing track: {bed_path}",
        verbosity,
        VERBOSITY_LVL[2],
    )
    try:
        if _is_bgzipped(bed_path):
            # bgzip-compressed but missing .tbi; no need to sort, just index
            print_verbosity(
                f"bgzip-compressed but missing .tbi index; building index only: {bed_path}",
                verbosity,
                VERBOSITY_LVL[3],
            )
            sorted_plain = bed_path
        else:
            # plain-text BED; must sort before compression
            print_verbosity(
                f"Plain-text BED; sorting before compression: {bed_path}",
                verbosity,
                VERBOSITY_LVL[3],
            )
            sorted_plain = _sort_bed(bed_path, debug)
        return _compress_and_index(sorted_plain, bed_path)
    except Exception as e:  # malformed BED, I/O error, indexing failure, ...
        exception_handler(
            CrispritzAnnotationError,
            f"Failed preparing annotation track {bed_path}",
            os.EX_DATAERR,
            debug,
            e,
        )


def _prepare_tracks(
    annotation_files: Sequence[str], threads: int, verbosity: int, debug: bool
) -> List[str]:
    """Prepare every annotation track for random-access querying.

    Calls :func:`_prepare_bed_track` on each file.  When more than one
    track is supplied and *threads* > 1, preparation is parallelised across
    tracks using a :class:`~concurrent.futures.ThreadPoolExecutor`.  The
    pool size is capped at ``min(threads, len(annotation_files))`` so no
    idle workers are created.

    The single-track case skips pool construction entirely to avoid the
    associated overhead.

    Parameters
    ----------
    annotation_files : Sequence[str]
        Paths to the raw BED annotation tracks.
    threads : int
        Maximum number of worker threads for parallel preparation.
    verbosity : int
        Controls progress message output via
        :func:`~crispritz_plus.verbosity.print_verbosity`.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    List[str]
        Paths to the bgzip-compressed, Tabix-indexed tracks in the same
        order as *annotation_files*.
    """
    workers = max(1, min(threads, len(annotation_files)))
    print_verbosity(
        f"Preparing {len(annotation_files)} track(s) using {workers} thread(s)",
        verbosity,
        VERBOSITY_LVL[2],
    )
    if workers == 1:  # avoid pool overhead for the common single-track case
        return [_prepare_bed_track(p, verbosity, debug) for p in annotation_files]
    prepared: List[Optional[str]] = [None] * len(annotation_files)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        index_to_future = {
            i: pool.submit(_prepare_bed_track, path, verbosity, debug)
            for i, path in enumerate(annotation_files)
        }
        for i, future in index_to_future.items():
            prepared[i] = future.result()
    return [p for p in prepared if p is not None]


# ==============================================================================
# Internal helpers - per-site annotation
# ==============================================================================


def _overlapping_features(tabix: pysam.TabixFile, offtarget: OffTarget) -> List[str]:
    """Return deduplicated BED feature labels overlapping an off-target site.

    Queries *tabix* for all records whose genomic interval intersects
    ``[offtarget.start, offtarget.end)``.  The fourth BED column (the
    feature label) is extracted from each overlapping record; duplicate
    labels within the same track are suppressed while preserving first-
    occurrence order.

    Parameters
    ----------
    tabix : pysam.TabixFile
        Open Tabix handle for one annotation track.
    offtarget : OffTarget
        Genomic footprint of the off-target site to query.

    Returns
    -------
    List[str]
        Ordered, deduplicated feature labels.  Returns ``[]`` when the
        chromosome is absent from the track index or no features overlap.
    """
    try:
        rows = tabix.fetch(offtarget.chrom, offtarget.start, offtarget.end)
    except ValueError:
        return []  # contig absent from this track's index -> no overlap
    labels: List[str] = []
    seen: set = set()
    for row in rows:
        fields = row.split("\t")
        label = fields[3]
        if label not in seen:
            seen.add(label)
            labels.append(label)
    return labels


def _genomic_footprint(spacer: str) -> int:
    """Return the number of non-gap bases in *spacer*.

    Counts characters that are **not** ``'-'`` to determine the number of
    reference bases consumed by the aligned spacer.  Gap characters (``'-'``)
    represent insertions in the guide relative to the reference and do not
    advance the genomic coordinate.

    Parameters
    ----------
    spacer : str
        Aligned spacer sequence, possibly containing ``'-'`` gap characters.

    Returns
    -------
    int
        Number of reference bases spanned by *spacer*.
    """
    return sum(base != "-" for base in spacer)


def _compute_target_footprint(pos1: int, spacer: str, strand: str) -> Tuple[int, int]:
    """Compute the 0-based half-open genomic interval of an off-target site.

    Converts the 1-based position from the targets TSV to a ``[start, end)``
    interval compatible with ``pysam.TabixFile.fetch``.

    Strand conventions
    ------------------
    * ``'+'`` (forward) - the stored position is the **leftmost** base of the
      spacer; the interval extends rightward by :func:`_genomic_footprint`
      bases.
    * ``'-'`` (reverse) - the stored position is the **rightmost** base of
      the spacer; the interval extends leftward by :func:`_genomic_footprint`
      bases.

    Parameters
    ----------
    pos1 : int
        1-based genomic position as stored in the targets TSV ``pos`` column.
    spacer : str
        Aligned spacer sequence (may contain ``'-'`` gap characters).
    strand : str
        Strand of the off-target site (``'+'`` or ``'-'``).

    Returns
    -------
    Tuple[int, int]
        ``(start, end)`` — a 0-based, half-open genomic interval.
    """
    pos0 = pos1 - 1  # 1-based -> 0-based half-open start
    if strand == "+":
        return pos0, pos0 + _genomic_footprint(spacer)
    return pos0 - _genomic_footprint(spacer), pos0


# ==============================================================================
# Internal helpers - streaming annotation loop
# ==============================================================================


def _update_header(
    header_line: str,
    fout: TextIO,
    column_names: List[str],
    targets_file: str,
    debug: bool,
) -> Dict[str, int]:
    """Validate the targets header, append annotation columns, and write it.

    Parses *header_line*, delegates schema validation to
    :func:`_validate_targets_header`, appends *column_names* to the right of
    the existing columns, and writes the extended header to *fout*.

    Parameters
    ----------
    header_line : str
        The raw (already stripped) first line of the targets TSV.
    fout : TextIO
        Output file object to which the extended header line is written.
    column_names : List[str]
        Annotation column names to append, one per BED track.
    targets_file : str
        Path to the targets file; used only to produce descriptive error
        messages.
    debug : bool
        When *True*, validation errors propagate with a full traceback.

    Returns
    -------
    Dict[str, int]
        Mapping of ``column_name -> 0-based index`` for all columns in the
        original (pre-annotation) header.

    Raises
    ------
    CrispritzAnnotationError
        If *header_line* is empty or the column schema is invalid.
    """
    if not header_line:
        exception_handler(
            CrispritzAnnotationError,
            f"Empty targets table: {targets_file}",
            os.EX_DATAERR,
            debug,
        )
    header_fields = header_line.split("\t")
    column_index = _validate_targets_header(header_fields, debug)
    fout.write("\t".join(header_fields + column_names) + "\n")
    return column_index


def _retrieve_offtarget_position(
    fields: List[str], idx: int, line_no: int, debug: bool
) -> int:
    """Extract and parse the 1-based genomic position from a TSV row.

    Parameters
    ----------
    fields : List[str]
        Tab-split fields of a single targets TSV data row.
    idx : int
        0-based column index of the ``pos`` field within *fields*.
    line_no : int
        1-based line number within the file (the header is line 1); used
        in error messages only.
    debug : bool
        When *True*, parse errors propagate with a full traceback.

    Returns
    -------
    int
        The parsed 1-based genomic position.

    Raises
    ------
    CrispritzAnnotationError
        If *idx* is out of range or the field value is not a valid integer.
    """
    try:
        return int(fields[idx])  # 1-based position (off-target convention)
    except (IndexError, ValueError) as e:
        line = "\t".join(fields)
        exception_handler(
            CrispritzAnnotationError,
            f"Malformed targets row at line {line_no}: {line!r}",
            os.EX_DATAERR,
            debug,
            e,
        )


def _read_offtarget(
    fields: List[str], column_index: Dict[str, int], line_no: int, debug: bool
) -> OffTarget:
    """Build an :class:`OffTarget` from the tab-split fields of a TSV row.

    Resolves chromosome, strand, and spacer from *fields* using
    *column_index*, then delegates coordinate conversion to
    :func:`_compute_target_footprint`.

    Parameters
    ----------
    fields : List[str]
        Tab-split fields of a single targets TSV data row.
    column_index : Dict[str, int]
        Mapping of ``column_name -> 0-based index`` as returned by
        :func:`_validate_targets_header`.
    line_no : int
        1-based line number within the file; used in error messages only.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    OffTarget
        Genomic footprint of the off-target site, ready for Tabix querying.

    Raises
    ------
    CrispritzAnnotationError
        If the ``pos`` field is missing or is not a valid integer.
    """
    pos1 = _retrieve_offtarget_position(fields, column_index["pos"], line_no, debug)
    chrom = fields[column_index["chrom"]]
    strand = fields[column_index["strand"]]
    spacer = fields[column_index["spacer"]]
    start, end = _compute_target_footprint(pos1, spacer, strand)
    return OffTarget(chrom=chrom, start=start, end=end, strand=strand)


def _annotate_stream(
    targets_file: str,
    fout: TextIO,
    tabix_handles: List[pysam.TabixFile],
    column_names: List[str],
    verbosity: int,
    debug: bool,
) -> int:
    """Stream-annotate the targets TSV, writing augmented rows to *fout*.

    Reads *targets_file* line-by-line, resolves each off-target site's
    0-based half-open genomic footprint via :func:`_read_offtarget`, queries
    every Tabix handle in *tabix_handles* for overlapping features, and emits
    one extended TSV row per site.

    Cell values follow these rules:

    * Multiple overlapping features within a single track are joined with
      a comma.
    * A track with no overlapping feature produces :data:`_NO_OVERLAP`
      (``"NA"``).

    Parameters
    ----------
    targets_file : str
        Path to the search-results TSV to annotate.
    fout : TextIO
        Writable output file object for the annotated TSV.
    tabix_handles : List[pysam.TabixFile]
        Open Tabix handles, one per annotation track, in the same order as
        *column_names*.
    column_names : List[str]
        Annotation column headers, one per track.
    verbosity : int
        Controls progress message output via
        :func:`~crispritz_plus.verbosity.print_verbosity`.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    int
        Total number of off-target rows annotated, excluding the header.
    """
    print_verbosity(
        f"Streaming annotation from {targets_file!r}",
        verbosity,
        VERBOSITY_LVL[3],
    )
    with open(targets_file, mode="r") as fin:
        header_line = fin.readline().strip()
        # Append the new headers to the right of the existing schema
        column_index = _update_header(
            header_line, fout, column_names, targets_file, debug
        )
        print_verbosity(
            f"Header validated: {len(column_index)} input column(s) + "
            f"{len(column_names)} annotation column(s)",
            verbosity,
            VERBOSITY_LVL[2],
        )
        annotated = 0
        for line_no, line in enumerate(fin, start=2):  # header was line 1
            if not (line := line.strip()):
                continue
            fields = line.split("\t")
            offtarget = _read_offtarget(fields, column_index, line_no, debug)
            cells = [
                (
                    ",".join(features)
                    if (features := _overlapping_features(tabix, offtarget))
                    else _NO_OVERLAP
                )
                for tabix in tabix_handles
            ]
            fout.write("\t".join(fields + cells) + "\n")
            annotated += 1
    print_verbosity(
        f"Annotated {annotated} off-target site(s)", verbosity, VERBOSITY_LVL[1]
    )
    return annotated


# ==============================================================================
# Public API
# ==============================================================================


def annotate_results(
    targets_file: str,
    annotation_files: List[str],
    outdir: str = os.getcwd(),
    annotation_names: Optional[List[str]] = None,
    threads: int = 1,
    verbosity: int = 1,
    debug: bool = False,
) -> None:
    """Annotate a CRISPRitz-plus search output with BED genomic feature tracks.

    This is the public entry point for the ``annotate-results`` CLI
    subcommand.  It orchestrates the full annotation pipeline in three phases:

    1. **Resolve column names** - derive output column headers from
       *annotation_names*, or generate ``annotation1`` … ``annotationN``
       defaults when *annotation_names* is ``None``.
    2. **Prepare tracks** - sort, bgzip-compress, and Tabix-index each BED
       file as needed; already-indexed files are passed through unchanged.
       Multiple tracks are prepared in parallel when *threads* > 1.
    3. **Stream-annotate** - iterate *targets_file* row-by-row, query each
       indexed track for overlapping genomic features, and write the augmented
       TSV to *outdir*.

    The output file is written to
    ``{outdir}/{stem_of_targets_file}.annotated.tsv``.

    Parameters
    ----------
    targets_file : str
        Path to the TSV produced by the ``search`` command.  The file must
        present the schema defined by :data:`SEARCH_OUTPUT_HEADER` in its
        first :data:`_CORE_COLUMNS` columns.
    annotation_files : List[str]
        One or more BED files providing genomic feature annotations.  Both
        plain-text and bgzip-compressed files are accepted; unindexed files
        are prepared automatically.
    outdir : str, optional
        Directory for the annotated output file.  Defaults to the current
        working directory at import time.
    annotation_names : Optional[List[str]], optional
        Custom output column names for each annotation track, in the same
        order as *annotation_files*.  When ``None``, columns are labelled
        ``annotation1`` … ``annotationN``.
    threads : int, optional
        Number of threads used when preparing annotation tracks in parallel.
        Defaults to ``1`` (single-threaded; no thread-pool overhead).
    verbosity : int, optional
        Controls the volume of progress output.  ``0`` suppresses all
        messages; higher values increase detail.  Defaults to ``1``.
    debug : bool, optional
        When *True*, exceptions propagate with a full traceback instead of
        being converted to a formatted user-facing error message.  Defaults
        to ``False``.

    Returns
    -------
    None
        Side-effect only: the annotated TSV is written to disk.

    Raises
    ------
    CrispritzAnnotationError
        On schema validation failure, I/O error, track preparation failure,
        or any other unexpected error during the annotation pipeline.

    Examples
    --------
    Annotate a search result with two BED tracks and explicit column names::

        from crispritz_plus.annotation import annotate_results

        annotate_results(
            targets_file="results/search_output.tsv",
            annotation_files=["tracks/genes.bed", "tracks/repeats.bed.gz"],
            outdir="results/",
            annotation_names=["genes", "repeats"],
            threads=4,
            verbosity=2,
        )
    """
    print_verbosity(
        f"annotate_results: targets={targets_file!r}, tracks={annotation_files}, "
        f"outdir={outdir!r}, threads={threads}",
        verbosity,
        VERBOSITY_LVL[3],
    )
    column_names = _resolve_annotation_names(annotation_files, annotation_names, debug)
    print_verbosity(
        f"Annotation column names: {column_names}",
        verbosity,
        VERBOSITY_LVL[2],
    )
    targets_ann_file = _create_targets_ann(targets_file, outdir)
    tabix_handles: List[pysam.TabixFile] = []
    try:
        print_verbosity(
            f"Preparing {len(annotation_files)} annotation track(s)",
            verbosity,
            VERBOSITY_LVL[1],
        )
        tabix_handles = [
            pysam.TabixFile(p)
            for p in _prepare_tracks(annotation_files, threads, verbosity, debug)
        ]
        print_verbosity(
            f"Annotating {targets_file} -> {targets_ann_file}",
            verbosity,
            VERBOSITY_LVL[1],
        )
        with open(targets_ann_file, "w") as fout:
            _annotate_stream(
                targets_file, fout, tabix_handles, column_names, verbosity, debug
            )
    except Exception as e:  # unexpected failure -> wrap consistently
        exception_handler(
            CrispritzAnnotationError,
            f"Failed annotating targets table {targets_file}",
            os.EX_SOFTWARE,
            debug,
            e,
        )
