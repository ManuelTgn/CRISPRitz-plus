"""
Plot-generation module for CRISPRitz-plus.

Produces two families of visualisations from an annotated search-results
TSV (the output of the ``annotate-results`` command):

Off-target sequence logos
    One PNG per guide RNA, rendered with LogoMaker_.  Each position in the
    logo shows the *deviation frequency* — how often off-target sites differ
    from the guide at that position.  Substitution mismatches are shown as
    the target base (A / C / G / T); DNA and RNA bulges are shown as the
    distinct symbols :data:`DNA_BULGE_SYMBOL` (``'D'``) and
    :data:`RNA_BULGE_SYMBOL` (``'R'``).

Annotation radar charts
    One polar chart per ``(annotation_column, mismatch_level)`` combination,
    for each guide and for the full cohort of guides.  Each apex of the radar
    represents one genomic feature label found in the annotation track; the
    radial value is the fraction of annotated off-targets that carry that
    label, relative to all labelled sites at the given mismatch ceiling.
    Features whose count falls below 10 % of the most frequent feature are
    suppressed to keep the chart legible.

Rendering note
--------------
``matplotlib`` is configured to use the ``Agg`` headless backend at import
time so that this module can be safely imported on servers and in CI
environments without a display.

Public API
----------
detect_annotation_columns
    Return the annotation-track column names present in a TSV header.
validate_annotated_tsv
    Validate that a TSV header satisfies the expected schema.
load_annotated_targets
    Read, validate, and type-cast an annotated search-results TSV.
build_offtarget_logo_matrix
    Build the per-position deviation-frequency matrix for one guide.
plot_offtarget_logo
    Render a LogoMaker sequence logo from a frequency matrix.
plot_annotation_radar
    Save one radar chart per ``(annotation_column, mismatch_level)``
    combination for a given subset of the data.
generate_all_plots
    Orchestrate all per-guide logos and annotation radars for a full report.

Module-level constants
----------------------
SEARCH_OUTPUT_HEADER : Tuple[str, ...]
    Fixed column names emitted by the ``search`` command.  Any column
    beyond this set is treated as an annotation track.
GAP_CHAR : str
    Gap / insertion character in aligned gRNA and spacer sequences
    (``'-'``).
DNA_BULGE_SYMBOL : str
    Logo-alphabet symbol for a DNA-bulge position (``'D'``).
RNA_BULGE_SYMBOL : str
    Logo-alphabet symbol for an RNA-bulge position (``'R'``).
LOGO_ALPHABET : Tuple[str, ...]
    Full six-character alphabet used in the per-position frequency matrix.
LOGO_COLOR_SCHEME : Dict[str, str]
    Hex colours assigned to each logo-alphabet symbol.
_ABSENT_TOKENS : frozenset
    Lowercase string tokens treated as "no annotation" / absent values
    by :func:`_is_present`.

.. _LogoMaker: https://logomaker.readthedocs.io
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib
import logomaker
import os

matplotlib.use("Agg")  # headless backend: report generation needs no display


from ..exception_handlers import exception_handler

from .crispritz_report_errors import CrispritzReportError


# ==============================================================================
# Data classes
# ==============================================================================


@dataclass
class FeaturesCounts:
    """Mapping of genomic feature labels to their off-target hit counts.

    Produced by :func:`_initialize_annotation_counts` for a single
    ``(annotation_column, mismatch_ceiling)`` slice of the data.  Features
    whose count falls below 10 % of the maximum observed count are excluded
    before construction to keep radar charts readable.

    Attributes
    ----------
    counts : Dict[str, int]
        Mapping of ``feature_label -> hit_count``.  Only features with a
        non-zero count above the 10 % threshold are present; ``'NA'``
        (no-overlap sentinel) is never included.
    """

    counts: Dict[str, int]


@dataclass
class MismatchesCounts:
    """Per-mismatch-level collection of :class:`FeaturesCounts`.

    Produced by :func:`_initialize_radar_counts`.  One entry exists for
    every mismatch level from ``0`` to ``max_mm`` (inclusive), enabling a
    full series of radar charts to be emitted from a single pass over the
    annotation data.

    Attributes
    ----------
    counts : Dict[int, FeaturesCounts]
        Mapping of ``mismatch_level -> FeaturesCounts`` where the
        ``FeaturesCounts`` for level *k* considers only off-target sites
        with ``mismatches <= k``.
    """

    counts: Dict[int, FeaturesCounts]


# ==============================================================================
# Module-level constants
# ==============================================================================

#: Fixed column names emitted by the ``search`` command and required in
#: every input file.  Any column beyond this set is treated as an
#: annotation track added by ``annotate-results``.
SEARCH_OUTPUT_HEADER: Tuple[str, ...] = (
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
)

#: Gap / insertion character used in aligned gRNA and spacer sequences.
#: A ``'-'`` in the gRNA field signals a DNA bulge; a ``'-'`` in the
#: spacer field signals an RNA bulge.
GAP_CHAR: str = "-"

#: Single-character symbol representing a **DNA-bulge** position in the
#: off-target logo alphabet.  A DNA bulge occurs when the target genome
#: carries an extra base that is absent from the guide RNA; the guide-body
#: cursor is *not* advanced at this position.
DNA_BULGE_SYMBOL: str = "D"

#: Single-character symbol representing an **RNA-bulge** position in the
#: off-target logo alphabet.  An RNA bulge occurs when the guide RNA carries
#: an extra base that is absent from the target genome; the guide-body
#: cursor *is* advanced at this position.
RNA_BULGE_SYMBOL: str = "R"

#: Full six-symbol alphabet used to build the per-position frequency matrix
#: passed to LogoMaker.  Contains the four canonical DNA bases plus the two
#: bulge symbols.
LOGO_ALPHABET: Tuple[str, ...] = (
    "A",
    "C",
    "G",
    "T",
    DNA_BULGE_SYMBOL,
    RNA_BULGE_SYMBOL,
)

#: Hex colour assigned to each symbol in the off-target sequence logo.
LOGO_COLOR_SCHEME: Dict[str, str] = {
    "A": "#2ca02c",
    "C": "#1f77b4",
    "G": "#f5b800",
    "T": "#d62728",
    DNA_BULGE_SYMBOL: "#9467bd",
    RNA_BULGE_SYMBOL: "#8c564b",
}

#: Lowercase string tokens considered equivalent to "absent / no annotation"
#: by :func:`_is_present`.  Any annotation value whose stripped, lower-cased
#: string representation is in this set is treated as a missing value.
_ABSENT_TOKENS = frozenset(
    {"", "n", "na", "n/a", ".", "nan", "none", "no", "0", "false", "-"}
)


# ==============================================================================
# Public helpers - TSV loading and validation
# ==============================================================================


def detect_annotation_columns(header: Sequence[str]) -> List[str]:
    """Return the annotation-track column names present in a TSV header.

    Any column name **not** in :data:`SEARCH_OUTPUT_HEADER` is considered an
    annotation column — i.e., a column added by the ``annotate-results``
    command.  Column order from *header* is preserved.

    Parameters
    ----------
    header : Sequence[str]
        Full list of column names from the TSV (e.g. ``list(df.columns)``).

    Returns
    -------
    List[str]
        Annotation column names in the order they appear in *header*.
        Returns an empty list when no annotation columns are present.
    """
    fixed = set(SEARCH_OUTPUT_HEADER)
    return [col for col in header if col not in fixed]


def validate_annotated_tsv(header: Sequence[str], debug: bool) -> List[str]:
    """Validate that a TSV header satisfies the annotated-search schema.

    Two conditions are checked in order:

    1. Every column in :data:`SEARCH_OUTPUT_HEADER` must be present (though
       not necessarily in position); missing columns indicate the file was
       not produced by the ``search`` command.
    2. At least one annotation column must be present; a TSV with only the
       fixed columns has not yet been annotated by ``annotate-results``.

    Parameters
    ----------
    header : Sequence[str]
        Full list of column names from the TSV.
    debug : bool
        When *True*, exceptions propagate with a full traceback instead of
        a formatted user-facing message.

    Returns
    -------
    List[str]
        Annotation column names (those beyond the fixed search schema),
        in header order.

    Raises
    ------
    CrispritzReportError
        If any required column is missing, or if no annotation columns are
        found.
    """
    missing = [col for col in SEARCH_OUTPUT_HEADER if col not in header]
    if missing:
        exception_handler(
            CrispritzReportError,
            "input is not a valid annotated search TSV; missing required "
            f"column(s): {', '.join(missing)}. Produce it with 'annotate-results'",
            os.EX_IOERR,
            debug,
        )
    annotation_columns = detect_annotation_columns(header)
    if not annotation_columns:
        exception_handler(
            CrispritzReportError,
            "input TSV contains only the fixed search columns and no annotation "
            "columns; it appears unannotated. Run 'annotate-results' first",
            os.EX_IOERR,
            debug,
        )
    return annotation_columns


def load_annotated_targets(tsv_path: str, debug: bool) -> pd.DataFrame:
    """Read, validate, and type-cast an annotated search-results TSV.

    All columns are initially read as strings (``dtype=str``) to avoid
    pandas' default NA-inference interfering with annotation sentinel values.
    After schema validation, the ``mismatches`` column is cast to ``int``.

    Parameters
    ----------
    tsv_path : str
        Path to the annotated TSV produced by ``annotate-results``.
    debug : bool
        When *True*, exceptions propagate with a full traceback instead of
        a formatted user-facing message.

    Returns
    -------
    pd.DataFrame
        Validated DataFrame with all columns as strings except
        ``'mismatches'``, which is cast to ``int``.

    Raises
    ------
    CrispritzReportError
        On I/O failure, parse error, empty file, schema validation failure,
        or a non-integer ``mismatches`` column.
    """
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, ValueError, pd.errors.ParserError) as e:
        exception_handler(
            CrispritzReportError,
            f"failed to read annotated TSV: {tsv_path}",
            os.EX_IOERR,
            debug,
            e,
        )
    if df.shape[1] == 0 or df.empty:
        exception_handler(
            CrispritzReportError,
            f"annotated TSV '{tsv_path}' contains no rows",
            os.EX_DATAERR,
            debug,
        )
    validate_annotated_tsv(list(df.columns), debug)
    try:
        df["mismatches"] = df["mismatches"].astype(int)
    except ValueError as e:
        exception_handler(
            CrispritzReportError,
            "column 'mismatches' is not integer-valued",
            os.EX_DATAERR,
            debug,
            e,
        )
    return df


# ==============================================================================
# Internal helpers - guide normalisation and name sanitisation
# ==============================================================================


def _guide_key(grna: str) -> str:
    """Return the gap-stripped normalised form of a guide RNA sequence.

    Removes all :data:`GAP_CHAR` (``'-'``) characters from *grna* to
    produce a canonical identifier that is invariant to bulge encoding.
    Used to group off-target rows by guide and to derive safe output
    filenames.

    Parameters
    ----------
    grna : str
        Raw gRNA sequence, possibly containing ``'-'`` gap characters.

    Returns
    -------
    str
        The input string with all ``'-'`` characters removed.
    """
    return grna.replace(GAP_CHAR, "")


def _safe_name(name: str) -> str:
    """Sanitise *name* for use as a component of a filesystem path.

    Replaces every non-alphanumeric character with an underscore.  Returns
    ``'guide'`` when the resulting string would be empty (e.g. *name* is
    all punctuation or whitespace).

    Parameters
    ----------
    name : str
        Arbitrary string to sanitise (typically a gap-stripped guide
        sequence).

    Returns
    -------
    str
        An alphanumeric-plus-underscore string safe for inclusion in a
        filename; never empty.
    """
    return "".join(ch if ch.isalnum() else "_" for ch in name) or "guide"


# ==============================================================================
# Internal helpers - frequency matrix construction
# ==============================================================================


def _compute_frequency_matrix(grnas: List[str], spacers: List[str]) -> pd.DataFrame:
    """Build a per-position deviation-frequency matrix for a set of off-targets.

    Each row corresponds to one position in the guide-body (gap characters
    excluded); each column corresponds to one symbol in :data:`LOGO_ALPHABET`.
    A cell value is the fraction of off-target sites that show the given
    deviation at that position.

    Position assignment rules
    -------------------------
    The aligned (grna, spacer) pairs are walked character by character using
    a *cursor* that tracks the current guide-body position:

    * **DNA bulge** (``grna[i] == GAP_CHAR``): the target has an extra base
      not present in the guide RNA.  The :data:`DNA_BULGE_SYMBOL` count at
      *cursor* is incremented, but *cursor* is **not** advanced (the guide
      body position is not consumed by the bulge).
    * **RNA bulge** (``spacer[i] == GAP_CHAR``): the guide RNA has an extra
      base not present in the target.  The :data:`RNA_BULGE_SYMBOL` count at
      *cursor* is incremented and *cursor* **is** advanced.
    * **Substitution mismatch**: neither character is a gap, neither is
      ``'N'``, and the guide and target bases differ.  The target base
      (``spacer[i].upper()``) count at *cursor* is incremented and *cursor*
      is advanced.
    * **Match or ambiguous base ('N')**: *cursor* is advanced; no count is
      added.

    Parameters
    ----------
    grnas : List[str]
        Aligned gRNA sequences (with gap characters) for all off-target rows.
    spacers : List[str]
        Aligned spacer / target sequences (with gap characters) for all
        off-target rows; must be the same length as *grnas*.

    Returns
    -------
    pd.DataFrame
        Shape ``(body_len, 6)`` DataFrame with columns
        :data:`LOGO_ALPHABET` and values in ``[0, 1]``, where
        ``body_len`` is the length of the longest gap-stripped guide.
    """
    body_len = max(len(_guide_key(g)) for g in grnas)
    sym_idx = {sym: i for i, sym in enumerate(LOGO_ALPHABET)}
    counts = np.zeros((body_len, len(LOGO_ALPHABET)), dtype=float)
    for grna, spacer in zip(grnas, spacers):
        cursor = 0
        for gch, tch in zip(grna, spacer):
            if cursor >= body_len:
                break
            if gch == GAP_CHAR:  # DNA bulge: target has an extra base
                counts[cursor, sym_idx[DNA_BULGE_SYMBOL]] += 1.0
                continue  # guide body position not consumed
            if tch == GAP_CHAR:  # RNA bulge: guide has an extra base
                counts[cursor, sym_idx[RNA_BULGE_SYMBOL]] += 1.0
                cursor += 1
                continue
            gup, tup = gch.upper(), tch.upper()
            if gup != "N" and tup != "N" and gup != tup and tup in sym_idx:
                counts[cursor, sym_idx[tup]] += 1.0  # substitution mismatch
            cursor += 1
    freq = counts / float(len(grnas))
    return pd.DataFrame(freq, columns=list(LOGO_ALPHABET))


def _strip_pam(matrix: pd.DataFrame, template: str) -> pd.DataFrame:
    """Remove PAM-encoded positions from a logo frequency matrix.

    Retains only the rows of *matrix* whose corresponding position in
    *template* is **not** ``'N'``.  The convention is that PAM positions are
    encoded as ``'N'`` in the guide template because they are not
    sequence-specific (e.g. the ``N`` in an ``NGG`` PAM); non-``'N'``
    positions are spacer positions that the logo should display.

    .. note::
       Fixed PAM bases (non-``'N'`` characters within the PAM region, such
       as the ``GG`` of an ``NGG`` PAM) are retained by this logic because
       they satisfy the keep condition.  Full PAM stripping would require
       explicit geometric metadata about the PAM length, which is not
       currently stored in the TSV.

    Parameters
    ----------
    matrix : pd.DataFrame
        Frequency matrix returned by :func:`_compute_frequency_matrix`.
    template : str
        Gap-stripped guide sequence used as a position mask.

    Returns
    -------
    pd.DataFrame
        Filtered matrix with a reset integer index.  If no positions pass
        the filter, the original *matrix* is returned unchanged.
    """
    n, _ = matrix.shape
    keep = [i for i in range(n) if i < len(template) and template[i].upper() != "N"]
    if keep:
        return pd.DataFrame(matrix.iloc[keep].reset_index(drop=True))
    return matrix


# ==============================================================================
# Public API - logo matrix and plot
# ==============================================================================


def build_offtarget_logo_matrix(
    df: pd.DataFrame,
    guide: Optional[str] = None,
    strip_pam: bool = True,
    debug: bool = False,
) -> Union[pd.Series, pd.DataFrame]:
    """Build the per-position deviation-frequency matrix for a guide's off-targets.

    Optionally filters *df* to a single guide, computes the frequency matrix
    via :func:`_compute_frequency_matrix`, and (when *strip_pam* is enabled
    and a *guide* is specified) removes PAM-encoded positions via
    :func:`_strip_pam`.  The matrix index is set to a 1-based
    ``'position'``-labelled :class:`~pandas.RangeIndex` for direct use by
    LogoMaker.

    Parameters
    ----------
    df : pd.DataFrame
        Full annotated targets DataFrame as returned by
        :func:`load_annotated_targets`.  Must contain ``'grna'`` and
        ``'spacer'`` columns.
    guide : Optional[str], optional
        Gap-stripped guide key to filter to.  When ``None`` the matrix is
        computed over all rows in *df*.
    strip_pam : bool, optional
        When ``True`` **and** *guide* is specified, remove positions encoded
        as ``'N'`` in the most-frequent guide template.  Defaults to
        ``True``.
    debug : bool, optional
        When ``True``, exceptions propagate with a full traceback.  Defaults
        to ``False``.

    Returns
    -------
    pd.DataFrame
        Shape ``(n_positions, 6)`` frequency matrix with columns
        :data:`LOGO_ALPHABET` and a 1-based ``'position'`` index, ready to
        pass to :func:`plot_offtarget_logo` or directly to
        ``logomaker.Logo``.

    Raises
    ------
    CrispritzReportError
        If *df* is empty (after optional guide filtering) so that no
        frequency matrix can be computed.
    """
    if guide is not None:
        df = pd.DataFrame(df[df["grna"].map(_guide_key) == guide])
    if df.empty:
        exception_handler(
            CrispritzReportError,
            "no off-targets available to build the off-target logo"
            + (f" for guide {guide!r}" if guide is not None else ""),
            os.EX_DATAERR,
            debug,
        )
    matrix = _compute_frequency_matrix(df["grna"].tolist(), df["spacer"].tolist())
    if strip_pam and guide is not None:
        template = df["grna"].map(_guide_key).mode().iloc[0]
        matrix = _strip_pam(matrix, template)
    matrix.index = pd.RangeIndex(start=1, stop=len(matrix) + 1, name="position")
    return matrix


def plot_offtarget_logo(
    matrix: pd.DataFrame,
    title: Optional[str] = None,
    ax: Optional[Axes] = None,
    debug: bool = False,
) -> Figure:
    """Render a LogoMaker sequence logo from a deviation-frequency matrix.

    Draws one coloured bar stack per guide position using :data:`LOGO_COLOR_SCHEME`.
    Only deviating symbols (mismatches or bulges) have non-zero height; perfect
    matches contribute nothing to the logo, so tall stacks indicate positions
    prone to off-target activity.

    A legend identifying the two bulge symbols is appended in the upper-right
    corner.

    Figure sizing
    -------------
    When *ax* is ``None`` a new :class:`~matplotlib.figure.Figure` is created
    whose width scales with the number of positions (``max(6.0, 0.45 x n)``
    inches) and whose height is fixed at 3 inches.  When *ax* is supplied
    the enclosing figure is reused.

    Parameters
    ----------
    matrix : pd.DataFrame
        Per-position frequency matrix as returned by
        :func:`build_offtarget_logo_matrix`.  Must have columns matching
        :data:`LOGO_ALPHABET` and at least one row.
    title : Optional[str], optional
        Axes title.  When ``None`` no title is set.
    ax : Optional[Axes], optional
        Pre-existing axes to draw into.  When ``None`` a new figure and axes
        are created.
    debug : bool, optional
        When ``True``, exceptions propagate with a full traceback.  Defaults
        to ``False``.

    Returns
    -------
    Figure
        The matplotlib figure containing the rendered logo.

    Raises
    ------
    CrispritzReportError
        If *matrix* is empty.
    """
    if matrix.empty or matrix.shape[0] == 0:
        exception_handler(
            CrispritzReportError,
            "cannot plot an empty off-target logo matrix",
            os.EX_DATAERR,
            debug,
        )
    if ax is None:
        width = max(6.0, 0.45 * len(matrix))
        fig, ax = plt.subplots(figsize=(width, 3.0))
    else:
        fig = ax.figure
    logo = logomaker.Logo(matrix, ax=ax, color_scheme=LOGO_COLOR_SCHEME)
    logo.style_spines(visible=False)
    logo.style_spines(spines=("left", "bottom"), visible=True)
    ax.set_xlabel("guide position")
    ax.set_ylabel("deviation frequency")
    if title:
        ax.set_title(title)
    # Legend distinguishing the two bulge symbols from substitutions.
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=LOGO_COLOR_SCHEME[DNA_BULGE_SYMBOL], label=f"DNA bulge ({DNA_BULGE_SYMBOL})"),  # type: ignore
        plt.Line2D([0], [0], marker="s", linestyle="", color=LOGO_COLOR_SCHEME[RNA_BULGE_SYMBOL], label=f"RNA bulge ({RNA_BULGE_SYMBOL})"),  # type: ignore
    ]
    ax.legend(handles=handles, loc="upper right", fontsize="small", frameon=False)
    fig.tight_layout()
    return fig


# ==============================================================================
# Internal helpers - radar chart data preparation
# ==============================================================================


def _initialize_annotation_counts(
    annotations: List[str], mismatches: List[int], max_mm: int
) -> FeaturesCounts:
    """Tally annotation-feature hits at a given mismatch ceiling.

    Counts how many off-target sites carry each annotation feature label
    (excluding the ``'NA'`` no-overlap sentinel) when the site's mismatch
    count does not exceed *max_mm*.  Features whose final count falls at or
    below 10 % of the maximum observed count are filtered out to prevent
    the radar chart from being cluttered by very rare features.

    Parameters
    ----------
    annotations : List[str]
        Annotation values from one track column, one entry per off-target
        row.  ``'NA'`` marks sites with no overlapping feature.
    mismatches : List[int]
        Mismatch counts for the corresponding off-target rows; must be the
        same length as *annotations*.
    max_mm : int
        Maximum mismatch count to include.  Sites with
        ``mismatches > max_mm`` are excluded from the tallies.

    Returns
    -------
    FeaturesCounts
        A :class:`FeaturesCounts` object containing only the features
        whose count is both positive and greater than 10 % of the
        maximum count.
    """
    counts = {f: 0 for f in set(annotations) if f != "NA"}
    for ann, mm in zip(annotations, mismatches):
        if ann != "NA" and mm <= max_mm:
            counts[ann] += 1
    thresh = max(list(counts.values())) * 0.1
    counts_ = {f: c for f, c in counts.items() if c > 0 and c > thresh}
    return FeaturesCounts(counts=counts_)


def _initialize_radar_counts(
    max_mm: Optional[int], mismatches: List[int], annotations: List[str]
) -> MismatchesCounts:
    """Build the full mismatch-stratified annotation counts for a radar series.

    Calls :func:`_initialize_annotation_counts` for every mismatch level
    from ``0`` up to *max_mm* (inclusive).  When *max_mm* is ``None``, the
    ceiling is inferred from the maximum value in *mismatches*.

    Parameters
    ----------
    max_mm : Optional[int]
        Maximum mismatch level to include.  When ``None``, the maximum
        value in *mismatches* is used (``0`` if *mismatches* is empty).
    mismatches : List[int]
        Mismatch counts, one per off-target row.
    annotations : List[str]
        Annotation values from one track column, one per off-target row;
        must be the same length as *mismatches*.

    Returns
    -------
    MismatchesCounts
        A :class:`MismatchesCounts` object mapping each mismatch level
        ``0 … max_mm`` to its :class:`FeaturesCounts`.
    """
    if max_mm is None:
        max_mm = max(mismatches, default=0)
    counts = {
        mm: _initialize_annotation_counts(annotations, mismatches, mm)
        for mm in range(max_mm + 1)
    }
    return MismatchesCounts(counts=counts)


# ==============================================================================
# Internal helpers - radar chart rendering and I/O
# ==============================================================================


def _plot_radar(
    features_counts: FeaturesCounts,
    mm: int,
    title: Optional[str] = None,
    color: str = "#1f77b4",
    fill_alpha: float = 0.25,
    label_pad: int = 14,
) -> Figure:
    """Render a polar radar chart for one mismatch level's annotation profile.

    Each apex of the radar represents one genomic feature label.  The radial
    value for each apex is the feature's share of all annotated off-target
    hits at the given mismatch ceiling (i.e. a proportion, not a raw count).
    Raw counts are appended to the apex labels in parentheses.

    Rendering details
    -----------------
    * The first apex is placed at the top (``θ_offset = π/2``).
    * The chart proceeds clockwise (``θ_direction = -1``).
    * Radial tick labels are suppressed; concentric rings are kept for
      reference.
    * The ``y`` axis upper limit is set to ``max(proportion) x 1.15`` to
      provide visual headroom.

    Parameters
    ----------
    features_counts : FeaturesCounts
        Annotation feature tallies for this mismatch level, as produced by
        :func:`_initialize_annotation_counts`.
    mm : int
        The mismatch ceiling represented by this chart; used in the default
        title string.
    title : Optional[str], optional
        Axes title override.  When ``None``, defaults to
        ``"Annotation profile (<mm> mismatches)"``.
    color : str, optional
        Hex colour for the radar line and fill.  Defaults to
        ``'#1f77b4'`` (Matplotlib's default blue).
    fill_alpha : float, optional
        Opacity of the filled radar polygon.  Defaults to ``0.25``.
    label_pad : int, optional
        Padding (points) between the chart edge and the apex labels.
        Defaults to ``14``.

    Returns
    -------
    Figure
        A 6 x 6 inch polar :class:`~matplotlib.figure.Figure`.
    """
    labels = list(features_counts.counts.keys())
    values = [float(v) for v in features_counts.counts.values()]
    labels_ = [f"{l} ({int(v)})" for l, v in zip(labels, values)]
    total = float(sum(values))
    values_pct = [v / total for v in values] if total > 0 else [0.0 for _ in values]
    n = len(labels)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angles_closed = np.concatenate([angles, angles[:1]])  # close the loop
    values_closed = values_pct + values_pct[:1]
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": "polar"})
    ax.set_theta_offset(np.pi / 2)  # first apex at top
    ax.set_theta_direction(-1)  # clockwise
    ax.set_thetagrids(np.degrees(angles), labels=labels_)
    ax.set_yticklabels([])  # suppress radial number labels (keep rings)
    ax.tick_params(axis="x", pad=label_pad)  # push apex labels outward
    ax.set_ylim(0, max(values_pct, default=1) * 1.15)
    ax.plot(angles_closed, values_closed, color=color, linewidth=2)
    ax.set_title(title or f"Annotation profile ({mm} mismatches)", pad=24)
    ax.fill(angles_closed, values_closed, color=color, alpha=fill_alpha)
    return fig


# ==============================================================================
# Public API - radar generation and figure I/O
# ==============================================================================


def plot_annotation_radar(
    targets_annotated: pd.DataFrame,
    annotation_columns: List[str],
    max_mm: Optional[int],
    title: str,
    dpi: int,
    outdir: str,
    prefix: str,
    outputs: List[str],
) -> List[str]:
    """Save one radar chart per ``(annotation_column, mismatch_level)`` pair.

    For each annotation track column in *annotation_columns*, computes a
    mismatch-stratified annotation profile via
    :func:`_initialize_radar_counts`, then renders and saves one PNG per
    mismatch level via :func:`_plot_radar` and :func:`_save_and_close`.

    Output filenames follow the pattern::

        {prefix}.{ann_col}.annotation_radar.{mm}.mismatches.png

    .. note::
       *outputs* is mutated **in place** (via ``list.append``) **and**
       returned, allowing both mutation-by-reference and fluent chaining
       from the caller.

    Parameters
    ----------
    targets_annotated : pd.DataFrame
        Annotated targets DataFrame (from :func:`load_annotated_targets`)
        containing ``'mismatches'`` (``int``) and all columns named in
        *annotation_columns*.
    annotation_columns : List[str]
        Names of the annotation track columns to iterate over.
    max_mm : Optional[int]
        Maximum mismatch ceiling for the radar series.  When ``None``, the
        maximum observed value in ``'mismatches'`` is used for each track.
    title : str
        Base title passed to :func:`_plot_radar` (the mismatch level is
        appended as ``"- {mm} mismatches"``).
    dpi : int
        Output resolution in dots-per-inch.
    outdir : str
        Directory to write PNG files into.
    prefix : str
        Filename prefix applied to every output figure.
    outputs : List[str]
        Accumulator list to which absolute output paths are appended.
        Modified in place.

    Returns
    -------
    List[str]
        The same list object as *outputs*, with any newly written paths
        appended.
    """
    mismatches: List[int] = targets_annotated["mismatches"].tolist()
    for ann_col in annotation_columns:
        mismatches_counts = _initialize_radar_counts(
            max_mm, mismatches, targets_annotated[ann_col].tolist()
        )
        for mm, features_counts in mismatches_counts.counts.items():
            title_ = f"{title} - {ann_col} - {mm} mismatches"
            fig = _plot_radar(features_counts, mm, title_)
            outputs.append(
                _save_and_close(
                    fig,
                    dpi,
                    outdir,
                    f"{prefix}.{ann_col}.annotation_radar.{mm}.mismatches.png",
                )
            )
    return outputs


def _save_and_close(fig: Figure, dpi: int, outdir: str, fname: str) -> str:
    """Save a matplotlib figure to disk and release its memory.

    Writes *fig* to ``{outdir}/{fname}`` at *dpi* resolution with
    ``bbox_inches='tight'`` so that axis labels are not clipped, then
    closes the figure to free the associated Agg renderer memory.

    Parameters
    ----------
    fig : Figure
        The matplotlib figure to save.
    dpi : int
        Output resolution in dots-per-inch.
    outdir : str
        Directory in which the file is written; must already exist.
    fname : str
        Filename (including extension) for the output file.

    Returns
    -------
    str
        Absolute path of the written file.
    """
    path = os.path.abspath(os.path.join(outdir, fname))
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def generate_all_plots(
    tsv_path: str,
    outdir: str,
    debug: bool,
    mm: Optional[int] = None,
    guides: Optional[Sequence[str]] = None,
    prefix: str = "report",
    dpi: int = 150,
) -> List[str]:
    """Orchestrate all plot generation for a complete ``generate-report`` run.

    Loads the annotated TSV, ensures *outdir* exists, then produces:

    1. **Per-guide off-target logos** - one PNG per guide key, written as::

           {prefix}.{safe_guide_name}.offtarget_logo.png

    2. **Overall annotation radar** - one radar per
       ``(annotation_column, mismatch_level)`` across all guides combined,
       written as::

           {prefix}.{ann_col}.annotation_radar.{mm}.mismatches.png

    3. **Per-guide annotation radars** - the same radar series repeated for
       each guide's off-target subset, written as::

           {prefix}.{safe_guide_name}.{ann_col}.annotation_radar.{mm}.mismatches.png

    Guide selection
    ---------------
    When *guides* is ``None``, guide keys are derived from the ``'grna'``
    column by stripping gap characters (via :func:`_guide_key`) and
    de-duplicating; the resulting set is sorted lexicographically.  When
    *guides* is provided explicitly, it is used as-is (no de-duplication or
    sorting is applied).

    Parameters
    ----------
    tsv_path : str
        Path to the annotated TSV produced by ``annotate-results``.
    outdir : str
        Directory for all output PNG files; created if absent.
    debug : bool
        When ``True``, exceptions propagate with a full traceback.
    mm : Optional[int], optional
        Maximum mismatch ceiling for radar charts.  When ``None``, the
        maximum observed mismatch count is used for each annotation column.
        Defaults to ``None``.
    guides : Optional[Sequence[str]], optional
        Explicit list of gap-stripped guide sequences to process.  When
        ``None``, all guide keys present in the TSV are processed.
        Defaults to ``None``.
    prefix : str, optional
        Filename prefix for every output PNG.  Defaults to ``'report'``.
    dpi : int, optional
        Output resolution in dots-per-inch.  Defaults to ``150``.

    Returns
    -------
    List[str]
        Absolute paths of all written PNG files, in generation order:
        logos first, then overall radars, then per-guide radars.

    Raises
    ------
    CrispritzReportError
        If the TSV is invalid, *guides* resolves to an empty list, or any
        individual plot fails.
    """
    df = load_annotated_targets(tsv_path, debug)
    annotation_columns = detect_annotation_columns(list(df.columns))
    os.makedirs(outdir, exist_ok=True)
    outputs: List[str] = []
    if guides is None:
        guide_keys = sorted({_guide_key(g) for g in df["grna"]})
    else:
        guide_keys = list(guides)
    if not guide_keys:
        exception_handler(
            CrispritzReportError,
            "no guides found in the annotated TSV",
            os.EX_DATAERR,
            debug,
        )
    # Per-guide off-target logos
    for gkey in guide_keys:
        matrix = build_offtarget_logo_matrix(df, guide=gkey, debug=debug)
        fig = plot_offtarget_logo(matrix, title=f"Off-target profile - {gkey}", debug=debug)  # type: ignore
        outputs.append(
            _save_and_close(
                fig, dpi, outdir, f"{prefix}.{_safe_name(gkey)}.offtarget_logo.png"
            )
        )
    # Overall annotation radar across all guides
    outputs = plot_annotation_radar(
        df,
        annotation_columns,
        mm,
        "Annotation profile - all guides",
        dpi,
        outdir,
        prefix,
        outputs,
    )
    for gkey in guide_keys:
        plot_annotation_radar(
            df,
            annotation_columns,
            mm,
            f"Annotation profile - {gkey}",
            dpi,
            outdir,
            f"{prefix}.{_safe_name(gkey)}",
            outputs,
        )
    return outputs
