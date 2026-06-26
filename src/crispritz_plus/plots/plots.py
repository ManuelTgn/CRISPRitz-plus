""" """

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


@dataclass
class FeaturesCounts:
    counts: Dict[str, int]


@dataclass
class MismatchesCounts:
    counts: Dict[int, FeaturesCounts]


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

GAP_CHAR: str = "-"

DNA_BULGE_SYMBOL: str = "D"

RNA_BULGE_SYMBOL: str = "R"

LOGO_ALPHABET: Tuple[str, ...] = (
    "A",
    "C",
    "G",
    "T",
    DNA_BULGE_SYMBOL,
    RNA_BULGE_SYMBOL,
)

LOGO_COLOR_SCHEME: Dict[str, str] = {
    "A": "#2ca02c",
    "C": "#1f77b4",
    "G": "#f5b800",
    "T": "#d62728",
    DNA_BULGE_SYMBOL: "#9467bd",
    RNA_BULGE_SYMBOL: "#8c564b",
}

_ABSENT_TOKENS = frozenset(
    {"", "n", "na", "n/a", ".", "nan", "none", "no", "0", "false", "-"}
)


def detect_annotation_columns(header: Sequence[str]) -> List[str]:
    fixed = set(SEARCH_OUTPUT_HEADER)
    return [col for col in header if col not in fixed]


def validate_annotated_tsv(header: Sequence[str], debug: bool) -> List[str]:
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


def _guide_key(grna: str) -> str:
    return grna.replace(GAP_CHAR, "")


def _is_present(value: object) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() not in _ABSENT_TOKENS


def _safe_name(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name) or "guide"


def _compute_frequency_matrix(grnas: List[str], spacers: List[str]) -> pd.DataFrame:
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
    n, _ = matrix.shape
    keep = [i for i in range(n) if i < len(template) and template[i].upper() != "N"]
    if keep:
        return pd.DataFrame(matrix.iloc[keep].reset_index(drop=True))
    return matrix


def build_offtarget_logo_matrix(
    df: pd.DataFrame,
    guide: Optional[str] = None,
    strip_pam: bool = True,
    debug: bool = False,
) -> Union[pd.Series, pd.DataFrame]:
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


def _initialize_annotation_counts(
    annotations: List[str], mismatches: List[int], max_mm: int
) -> FeaturesCounts:
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
    if max_mm is None:
        max_mm = max(mismatches, default=0)
    counts = {
        mm: _initialize_annotation_counts(annotations, mismatches, mm)
        for mm in range(max_mm + 1)
    }
    return MismatchesCounts(counts=counts)


def _plot_radar(
    features_counts: FeaturesCounts,
    mm: int,
    title: Optional[str] = None,
    color: str = "#1f77b4",
    fill_alpha: float = 0.25,
    label_pad: int = 14,
) -> Figure:
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
    mismatches: List[int] = targets_annotated["mismatches"].tolist()
    for ann_col in annotation_columns:
        mismatches_counts = _initialize_radar_counts(
            max_mm, mismatches, targets_annotated[ann_col].tolist()
        )
        for mm, features_counts in mismatches_counts.counts.items():
            title_ = f"{title} - {mm} mismatches"
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
