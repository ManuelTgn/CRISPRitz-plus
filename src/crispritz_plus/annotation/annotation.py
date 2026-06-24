"""
Off-target annotation against genomic BED interval tracks.

This module is the modern, importable replacement for the legacy
``sourceCode/Python_Scripts/Annotator/annotator.py`` script. It consumes the
tab-separated targets table produced by the ``search`` command and, for each
supplied BED annotation track, appends one column reporting which features (if
any) overlap every off-target site.

Design notes
------------
* Overlap is resolved with :class:`pysam.TabixFile` and its ``.fetch()``
  method (block-gzip + tabix random access), so the genome-wide targets table
  is streamed once and every row triggers an O(log n) coordinate lookup per
  track rather than building an in-memory interval tree over the whole
  annotation (the legacy approach).
* Both uncompressed ``.bed`` and block-gzipped ``.bed.gz`` inputs are
  accepted. Any track that is not already a tabix-indexed ``.bed.gz`` is
  normalised (coordinate-sorted, bgzipped, indexed) into a private temporary
  directory; the user's original files are never modified, and every
  temporary artefact is removed when annotation finishes or errors out.
* The targets table is validated against the canonical ``search`` output
  schema before any work begins, so a wrong input fails fast with a clear,
  user-facing message.
"""



from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

import contextlib
import gzip
import os
import pysam
import shutil
import tempfile


from ..crispritz_errors import CrispritzError
from ..exception_handlers import exception_handler
from ..verbosity import VERBOSITY_LVL, print_verbosity
from .crispritz_annotation_error import CrispritzAnnotationError


# =============================================================================
# Canonical search-output schema
# =============================================================================

#: Column layout of the final targets table emitted by the ``search`` command.
#: This mirrors the C++ ``ScoredTsvFormatter`` / ``scores.shard_scoring``
#: contract and the merged-table schema in ``result_merger.hpp`` — NOT the
#: separately tracked, known-divergent ``TSV_HEADER`` constant in
#: ``offtarget.py``. The trailing ``cfd_score`` column is always present in a
#: merged search table, but it is treated as optional here so that unscored
#: intermediate tables still validate.
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

#: Columns required regardless of whether CFD scoring has run.
_CORE_COLUMNS: int = 9

#: Cell value written when a site overlaps no feature in a given track.
_NO_OVERLAP: str = "NA"

#: Suffix of the file produced by this command.
_OUTPUT_SUFFIX: str = ".annotated.tsv"

#: BED comment / metadata line prefixes to skip while normalising a track.
_BED_SKIP_PREFIXES: Tuple[str, ...] = ("#", "track", "browser")


# =============================================================================
# Header / argument validation
# =============================================================================


def _validate_targets_header(
    header_fields: Sequence[str], debug: bool
) -> Dict[str, int]:
    """Validate that *header_fields* is a ``search`` targets-table header.

    Parameters
    ----------
    header_fields:
        Tab-split fields of the input file's first line.
    debug:
        When *True*, validation failures propagate with a full traceback.

    Returns
    -------
    Dict[str, int]
        Mapping from canonical column name to its 0-based index, so downstream
        code can resolve fields by name rather than by a hard-coded position.
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
    annotation_files:
        The BED tracks supplied on the command line.
    annotation_names:
        Optional explicit names, one per track. When ``None`` the columns
        default to ``annotation1`` … ``annotationN``.
    debug:
        When *True*, a length mismatch propagates with a full traceback.

    Returns
    -------
    List[str]
        One column name per annotation track, in input order.
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


# =============================================================================
# BED track preparation (sort / bgzip / tabix-index)
# =============================================================================


def _is_bgzipped(path: str) -> bool:
    """Return whether *path* names a (b)gzip-compressed file by extension."""
    return path.endswith(".gz")


def _sort_bed(bed_path: str, tmpdir: str, debug: bool) -> str:
    records: List[Tuple[str, int, str]] = []
    try:
        opener = gzip.open if _is_bgzipped(bed_path) else open
        with opener(bed_path, mode="rt") as fin:
            for line in fin:
                line = line.rstrip("\n")
                if not line or line.startswith(_BED_SKIP_PREFIXES):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:  # assume annotation on 4th column
                    continue  # not a valid BED interval record
                records.append((parts[0], int(parts[1]), line))
        records.sort(key=lambda rec: (rec[0], rec[1], rec[2]))
        sorted_plain = os.path.join(tmpdir, f"{os.path.basename(bed_path)}.sorted.bed")
        with open(sorted_plain, "w") as fout:
            fout.writelines(f"{line}\n" for _, _, line in records)
    except (IOError, Exception) as e:
        exception_handler(CrispritzAnnotationError, f"Failed sorting BED {bed_path}", os.EX_IOERR, debug, e)
    return sorted_plain


def _prepare_bed_track(bed_path: str, tmpdir: str, debug: bool) -> str:
    """Return a path to a tabix-indexed, bgzipped, coordinate-sorted BED.

    If *bed_path* is already a bgzipped BED accompanied by a sidecar ``.tbi``
    index it is used in place. Otherwise the track is normalised into *tmpdir*:
    its records are coordinate-sorted (tabix indexing requires position-sorted
    input), block-gzipped, and tabix-indexed. The user's original file is never
    modified.

    Parameters
    ----------
    bed_path:
        Path to the user-supplied ``.bed`` or ``.bed.gz`` track.
    tmpdir:
        Private temporary directory that receives any generated artefacts.
    debug:
        When *True*, preparation failures propagate with a full traceback.

    Returns
    -------
    str
        Path to a tabix-indexed ``.bed.gz`` ready for ``TabixFile.fetch``.
    """
    if _is_bgzipped(bed_path) and os.path.isfile(f"{bed_path}.tbi"):
        return bed_path  # already random-access ready; nothing to do
    try:
        # Read (transparently decompressing a .gz) and coordinate-sort. gzip can
        # read a block-gzip stream, so the same path covers an indexless .bed.gz.
        sorted_plain = _sort_bed(bed_path, tmpdir, debug)
        gz_path = f"{sorted_plain}.gz"
        pysam.tabix_compress(sorted_plain, gz_path, force=True)
        pysam.tabix_index(gz_path, preset="bed", force=True)  # writes gz_path.tbi
        return gz_path
    except Exception as e:  # malformed BED, I/O error, indexing failure, ...
        exception_handler(
            CrispritzAnnotationError,
            f"Failed preparing annotation track {bed_path}",
            os.EX_DATAERR,
            debug,
            e,
        )


def _prepare_tracks(
    annotation_files: Sequence[str], tmpdir: str, threads: int, debug: bool
) -> List[str]:
    """Prepare every BED track, optionally in parallel, preserving order.

    Order preservation is essential: the prepared-track order maps one-to-one
    onto the appended column order and therefore onto the resolved column names.

    Parameters
    ----------
    annotation_files:
        BED tracks to prepare.
    tmpdir:
        Private temporary directory for generated artefacts.
    threads:
        Maximum worker threads (track preparation is independent per file and
        spends its time in GIL-releasing pysam C calls).
    debug:
        When *True*, failures propagate with a full traceback.

    Returns
    -------
    List[str]
        Tabix-ready ``.bed.gz`` paths, in the same order as *annotation_files*.
    """
    workers = max(1, min(threads, len(annotation_files)))
    if workers == 1:  # avoid pool overhead for the common single-track case
        return [_prepare_bed_track(p, tmpdir, debug) for p in annotation_files]
    prepared: List[Optional[str]] = [None] * len(annotation_files)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        index_to_future = {
            i: pool.submit(_prepare_bed_track, path, tmpdir, debug)
            for i, path in enumerate(annotation_files)
        }
        for i, future in index_to_future.items():
            prepared[i] = future.result()
    return [p for p in prepared if p is not None]


# =============================================================================
# Overlap resolution
# =============================================================================


def _overlapping_features(
    tabix: pysam.TabixFile, chrom: str, start: int, end: int
) -> List[str]:
    """Return the de-duplicated feature labels overlapping ``[start, end)``.

    A feature label is the BED name column (field 4) when present, otherwise
    the ``chrom:start-end`` region string. Coordinates are 0-based half-open,
    matching both BED and :meth:`pysam.TabixFile.fetch`.

    Parameters
    ----------
    tabix:
        An open tabix handle for one prepared BED track.
    chrom:
        Contig of the off-target site.
    start:
        0-based inclusive start of the query window.
    end:
        0-based exclusive end of the query window.

    Returns
    -------
    List[str]
        Overlapping feature labels in first-seen order (empty when none).
    """
    try:
        rows = tabix.fetch(chrom, start, end)
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
    """Return the reference span of an aligned spacer (gap characters excluded).

    Bulges introduce ``'-'`` gap characters into the aligned spacer; only the
    non-gap bases consume reference coordinates, so they alone define the
    genomic window queried against each annotation track.
    """
    return sum(base != "-" for base in spacer)


def _compute_target_footprint(pos1: int, spacer: str, strand: str) -> Tuple[int, int]:
    pos0 = pos1 - 1  # 1-based -> 0-based half-open start
    if strand == "+":
        return pos0, pos0 + _genomic_footprint(spacer) 
    return pos0 - _genomic_footprint(spacer), pos0


# =============================================================================
# Streaming annotation
# =============================================================================

def _update_header(header_line: str, column_names: List[str], targets_file: str, debug: bool) -> Tuple[Dict[str, int], List[str]]:
    if not header_line:
        exception_handler(CrispritzAnnotationError, f"Empty targets table: {targets_file}", os.EX_DATAERR, debug)
    header_fields = header_line.split("\t")
    column_index = _validate_targets_header(header_fields, debug)
    # Append the new headers to the right of the existing schema.
    return column_index, header_fields + column_names


def _annotate_stream(
    targets_file: str,
    fout: TextIO,
    tabix_handles: List[pysam.TabixFile],
    column_names: List[str],
    verbosity: int,
    debug: bool,
) -> int:
    """Stream the targets table, appending one annotation column per track.

    The input is read line-by-line and the output is written line-by-line, so
    peak memory is bounded by a single row regardless of table size.

    Parameters
    ----------
    targets_file:
        Path to the ``search`` targets table.
    fout:
        Open, writable text handle for the annotated output.
    tabix_handles:
        Open tabix handles, one per track, in column order.
    column_names:
        Output header names for the appended columns, in the same order.
    verbosity:
        Verbosity level (see :data:`crispritz_plus.verbosity.VERBOSITY_LVL`).
    debug:
        When *True*, failures propagate with a full traceback.

    Returns
    -------
    int
        Number of off-target rows annotated.
    """
    with open(targets_file, "r") as fin:
        header_line = fin.readline().rstrip("\n")
        column_index, header_fields = _update_header(header_line, column_names, targets_file, debug)
        # Append the new headers to the right of the existing schema.
        fout.write("\t".join(header_fields) + "\n")
        annotated = 0
        for line_no, raw in enumerate(fin, start=2):  # header was line 1
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            try:
                pos1 = int(fields[column_index["pos"]])  # 1-based position (OffTarget convention)
            except (IndexError, ValueError) as e:
                exception_handler(
                    CrispritzAnnotationError,
                    f"Malformed targets row at line {line_no}: {line!r}",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
            chrom = fields[column_index["chrom"]]
            start, end = _compute_target_footprint(pos1, fields[column_index["spacer"]], fields[column_index["strand"]])
            cells = [
                ",".join(features) if (features := _overlapping_features(
                    tabix, chrom, start, end)) else _NO_OVERLAP
                for tabix in tabix_handles
            ]
            fout.write("\t".join(fields + cells) + "\n")
            annotated += 1
    print_verbosity(
        f"Annotated {annotated} off-target site(s)", verbosity, VERBOSITY_LVL[1]
    )
    return annotated


# =============================================================================
# Public entry point
# =============================================================================


def annotate_results(
    targets_file: str,
    annotation_files: List[str],
    outdir: str = os.getcwd(),
    annotation_names: Optional[List[str]] = None,
    threads: int = 1,
    verbosity: int = 1,
    debug: bool = False,
) -> str:
    """Annotate a search targets table with one column per BED track.

    For every annotation track a column is appended to the right of the input
    targets table, listing the features that overlap each off-target site (or
    ``'.'`` when none overlap). All BED handling is performed against a private
    temporary directory that is removed unconditionally on completion or error.

    Parameters
    ----------
    targets_file:
        Path to the TSV targets table produced by the ``search`` command.
    annotation_files:
        One or more BED tracks (``.bed`` or block-gzipped ``.bed.gz``).
    outdir:
        Directory the annotated table is written to. The output file is named
        ``<targets-stem>.annotated.tsv``.
    annotation_names:
        Optional output-column names, one per track and in the same order. When
        omitted, columns default to ``annotation1`` … ``annotationN``.
    threads:
        Maximum worker threads used to prepare (sort/bgzip/index) BED tracks.
    verbosity:
        Verbosity level (see :data:`crispritz_plus.verbosity.VERBOSITY_LVL`).
    debug:
        When *True*, exceptions propagate with full tracebacks.

    Returns
    -------
    str
        Path to the annotated targets table.

    Raises
    ------
    CrispritzAnnotationError
        If the input is not a valid search targets table, a track cannot be
        prepared, or the annotation pass fails (full traceback only in debug
        mode; otherwise a clean message is printed and the process exits).
    """
    column_names = _resolve_annotation_names(annotation_files, annotation_names, debug)
    stem = os.path.basename(targets_file).rsplit(".", 1)[0]
    output_file = os.path.join(os.path.abspath(outdir), f"{stem}{_OUTPUT_SUFFIX}")
    tmpdir = tempfile.mkdtemp(prefix="crispritz_annotate_")
    tabix_handles: List[pysam.TabixFile] = []
    try:
        print_verbosity(f"Preparing {len(annotation_files)} annotation track(s)", verbosity, VERBOSITY_LVL[2])
        tabix_handles = [pysam.TabixFile(p) for p in _prepare_tracks(annotation_files, tmpdir, threads, debug)]
        print_verbosity(f"Annotating {targets_file} -> {output_file}", verbosity, VERBOSITY_LVL[1])
        with open(output_file, "w") as fout:
            _annotate_stream(
                targets_file,
                fout,
                tabix_handles,
                column_names,
                verbosity,
                debug,
            )
    except CrispritzError:
        raise  # already a clean, domain-specific error (debug mode)
    except Exception as e:  # unexpected failure -> wrap consistently
        exception_handler(
            CrispritzAnnotationError,
            f"Failed annotating targets table {targets_file}",
            os.EX_SOFTWARE,
            debug,
            e,
        )
    finally:
        # Strict cleanup: close handles and delete every temporary artefact,
        # whether annotation completed or raised.
        for handle in tabix_handles:
            with contextlib.suppress(Exception):
                handle.close()
        shutil.rmtree(tmpdir, ignore_errors=True)
    return output_file