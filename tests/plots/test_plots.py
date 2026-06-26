"""
tests/plots/test_plots.py
=========================
Unit tests for ``crispritz_plus.plots.plots``.

Scope
-----
Covers pure-logic helpers and public API functions.  Rendering fidelity
(visual correctness of logos/radar charts) is out of scope; tests that
touch matplotlib only verify that a ``Figure`` is returned and files are
written — not their pixel content.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt

from crispritz_plus.plots.plots import (
    # Constants
    SEARCH_OUTPUT_HEADER,
    GAP_CHAR,
    DNA_BULGE_SYMBOL,
    RNA_BULGE_SYMBOL,
    LOGO_ALPHABET,
    LOGO_COLOR_SCHEME,
    _ABSENT_TOKENS,
    # Data classes
    FeaturesCounts,
    MismatchesCounts,
    # Public helpers
    detect_annotation_columns,
    validate_annotated_tsv,
    load_annotated_targets,
    build_offtarget_logo_matrix,
    plot_offtarget_logo,
    plot_annotation_radar,
    # Private helpers
    _guide_key,
    _safe_name,
    _compute_frequency_matrix,
    _strip_pam,
    _initialize_annotation_counts,
    _initialize_radar_counts,
    _save_and_close,
)
from crispritz_plus.plots.crispritz_report_errors import CrispritzReportError


# ===========================================================================
# Module-level constants
# ===========================================================================

class TestConstants:
    def test_logo_alphabet_has_six_symbols(self):
        assert len(LOGO_ALPHABET) == 6

    def test_logo_alphabet_contains_bulge_symbols(self):
        assert DNA_BULGE_SYMBOL in LOGO_ALPHABET
        assert RNA_BULGE_SYMBOL in LOGO_ALPHABET

    def test_color_scheme_keys_match_alphabet(self):
        assert set(LOGO_COLOR_SCHEME.keys()) == set(LOGO_ALPHABET)

    def test_absent_tokens_is_frozenset(self):
        assert isinstance(_ABSENT_TOKENS, frozenset)

    def test_search_output_header_length(self):
        assert len(SEARCH_OUTPUT_HEADER) == 10


# ===========================================================================
# detect_annotation_columns
# ===========================================================================

class TestDetectAnnotationColumns:
    def test_no_extras_returns_empty(self):
        assert detect_annotation_columns(list(SEARCH_OUTPUT_HEADER)) == []

    def test_extra_column_detected(self):
        header = list(SEARCH_OUTPUT_HEADER) + ["gene_track"]
        assert detect_annotation_columns(header) == ["gene_track"]

    def test_multiple_extras_in_order(self):
        header = list(SEARCH_OUTPUT_HEADER) + ["genes", "repeats"]
        assert detect_annotation_columns(header) == ["genes", "repeats"]

    def test_empty_header_returns_empty(self):
        assert detect_annotation_columns([]) == []


# ===========================================================================
# validate_annotated_tsv
# ===========================================================================

class TestValidateAnnotatedTsv:
    def test_valid_header_returns_annotation_cols(self):
        header = list(SEARCH_OUTPUT_HEADER) + ["gene_track"]
        result = validate_annotated_tsv(header, debug=False)
        assert result == ["gene_track"]

    def test_missing_required_column_raises(self):
        # Remove "chrom" from header
        header = [c for c in SEARCH_OUTPUT_HEADER if c != "chrom"] + ["gene_track"]
        with pytest.raises(CrispritzReportError):
            validate_annotated_tsv(header, debug=False)

    def test_no_annotation_columns_raises(self):
        with pytest.raises(CrispritzReportError):
            validate_annotated_tsv(list(SEARCH_OUTPUT_HEADER), debug=False)


# ===========================================================================
# load_annotated_targets
# ===========================================================================

class TestLoadAnnotatedTargets:
    def test_valid_tsv_returns_dataframe(self, annotated_tsv):
        df = load_annotated_targets(annotated_tsv, debug=False)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty

    def test_mismatches_column_cast_to_int(self, annotated_tsv):
        df = load_annotated_targets(annotated_tsv, debug=False)
        assert df["mismatches"].dtype == int

    def test_missing_file_raises(self):
        with pytest.raises(CrispritzReportError):
            load_annotated_targets("/nonexistent/path.tsv", debug=False)

    def test_non_integer_mismatches_raises(self, tmp_path):
        header = "\t".join(list(SEARCH_OUTPUT_HEADER) + ["genes"])
        row    = "chr1\t100\t+\tACGT\tACGT\tNOT_INT\tX\t0\t0\t1.0\tprom"
        p = tmp_path / "bad.tsv"
        p.write_text(f"{header}\n{row}\n")
        with pytest.raises(CrispritzReportError):
            load_annotated_targets(str(p), debug=False)


# ===========================================================================
# _guide_key
# ===========================================================================

class TestGuideKey:
    def test_no_gaps_unchanged(self):
        assert _guide_key("ACGT") == "ACGT"

    def test_gaps_stripped(self):
        assert _guide_key("AC-GT") == "ACGT"

    def test_all_gaps_returns_empty(self):
        assert _guide_key("----") == ""


# ===========================================================================
# _safe_name
# ===========================================================================

class TestSafeName:
    def test_alphanumeric_unchanged(self):
        assert _safe_name("ACGT123") == "ACGT123"

    def test_special_chars_replaced(self):
        assert _safe_name("AC GT!") == "AC_GT_"

    def test_empty_result_becomes_guide(self):
        assert _safe_name("!@#") == "___"   # replaced; non-empty result

    def test_all_punctuation_returns_underscores_not_empty(self):
        result = _safe_name("--")
        assert result  # never empty


# ===========================================================================
# _compute_frequency_matrix
# ===========================================================================

class TestComputeFrequencyMatrix:
    def test_output_shape(self):
        grnas   = ["ACGT"]
        spacers = ["ACGT"]
        m = _compute_frequency_matrix(grnas, spacers)
        assert m.shape == (4, len(LOGO_ALPHABET))

    def test_perfect_match_all_zeros(self):
        m = _compute_frequency_matrix(["ACGT"], ["ACGT"])
        assert (m.values == 0).all()

    def test_substitution_mismatch_recorded(self):
        # pos 1: C→T mismatch
        m = _compute_frequency_matrix(["ACGT"], ["ATGT"])
        assert m.loc[1, "T"] == pytest.approx(1.0)

    def test_dna_bulge_recorded(self):
        # grna has gap at pos 0 → DNA bulge at cursor 0
        m = _compute_frequency_matrix(["-CGT"], ["ACGT"])
        assert m.loc[0, DNA_BULGE_SYMBOL] == pytest.approx(1.0)

    def test_rna_bulge_recorded(self):
        # spacer has gap at pos 0 → RNA bulge at cursor 0
        m = _compute_frequency_matrix(["ACGT"], ["-CGT"])
        assert m.loc[0, RNA_BULGE_SYMBOL] == pytest.approx(1.0)

    def test_frequency_normalised_by_n_rows(self):
        # Two identical mismatches → frequency stays 1.0 (normalised by 2)
        m = _compute_frequency_matrix(["ACGT", "ACGT"], ["ATGT", "ATGT"])
        assert m.loc[1, "T"] == pytest.approx(1.0)


# ===========================================================================
# _strip_pam
# ===========================================================================

class TestStripPam:
    def _matrix(self, n_rows: int) -> pd.DataFrame:
        return pd.DataFrame(
            np.zeros((n_rows, len(LOGO_ALPHABET))), columns=list(LOGO_ALPHABET)
        )

    def test_n_positions_removed(self):
        m = self._matrix(4)
        result = _strip_pam(m, "NACT")  # pos 0 is N → 3 rows kept
        assert len(result) == 3

    def test_non_n_positions_kept(self):
        m = self._matrix(3)
        result = _strip_pam(m, "ACT")  # no N → all kept
        assert len(result) == 3

    def test_all_n_returns_original(self):
        m = self._matrix(3)
        result = _strip_pam(m, "NNN")
        assert result is m  # unchanged original returned

    def test_index_reset_after_strip(self):
        m = self._matrix(4)
        result = _strip_pam(m, "NACT")
        assert list(result.index) == [0, 1, 2]


# ===========================================================================
# build_offtarget_logo_matrix
# ===========================================================================

class TestBuildOfftargetLogoMatrix:
    def test_returns_dataframe_with_alphabet_columns(self, minimal_df):
        m = build_offtarget_logo_matrix(minimal_df, guide=None, strip_pam=False)
        assert list(m.columns) == list(LOGO_ALPHABET)

    def test_guide_filter_applied(self, minimal_df):
        m = build_offtarget_logo_matrix(minimal_df, guide="ACGTACGT", strip_pam=False)
        assert isinstance(m, pd.DataFrame)
        assert not m.empty

    def test_unknown_guide_raises(self, minimal_df):
        with pytest.raises(CrispritzReportError):
            build_offtarget_logo_matrix(minimal_df, guide="ZZZZZZZZ", strip_pam=False)

    def test_empty_df_raises(self):
        empty = pd.DataFrame(columns=["grna", "spacer"])
        with pytest.raises(CrispritzReportError):
            build_offtarget_logo_matrix(empty, guide=None, strip_pam=False)


# ===========================================================================
# _initialize_annotation_counts
# ===========================================================================

class TestInitializeAnnotationCounts:
    def test_na_values_excluded(self):
        ann = ["promoter", "NA", "exon"]
        mm  = [0, 0, 0]
        fc  = _initialize_annotation_counts(ann, mm, max_mm=1)
        assert "NA" not in fc.counts

    def test_sites_above_max_mm_excluded(self):
        ann = ["promoter", "exon"]
        mm  = [0, 2]
        fc  = _initialize_annotation_counts(ann, mm, max_mm=1)
        # exon is at mm=2 which is > max_mm=1 → excluded
        assert fc.counts.get("exon", 0) == 0

    def test_rare_features_below_threshold_excluded(self):
        # promoter appears 10×, rare appears 1× → rare below 10% threshold
        ann = ["promoter"] * 10 + ["rare"]
        mm  = [0] * 11
        fc  = _initialize_annotation_counts(ann, mm, max_mm=0)
        assert "rare" not in fc.counts
        assert fc.counts["promoter"] == 10

    def test_returns_features_counts_instance(self):
        fc = _initialize_annotation_counts(["promoter"], [0], max_mm=0)
        assert isinstance(fc, FeaturesCounts)


# ===========================================================================
# _initialize_radar_counts
# ===========================================================================

class TestInitializeRadarCounts:
    def test_keys_range_zero_to_max_mm(self):
        mc = _initialize_radar_counts(3, [0, 1, 2, 3], ["A", "B", "C", "D"])
        assert set(mc.counts.keys()) == {0, 1, 2, 3}

    def test_max_mm_none_inferred_from_mismatches(self):
        mc = _initialize_radar_counts(None, [0, 1, 2], ["A", "B", "C"])
        assert set(mc.counts.keys()) == {0, 1, 2}

    def test_returns_mismatches_counts_instance(self):
        mc = _initialize_radar_counts(1, [0, 1], ["A", "B"])
        assert isinstance(mc, MismatchesCounts)


# ===========================================================================
# _save_and_close
# ===========================================================================

class TestSaveAndClose:
    def test_file_created_on_disk(self, tmp_path):
        fig, _ = plt.subplots()
        _save_and_close(fig, 72, str(tmp_path), "out.png")
        assert (tmp_path / "out.png").is_file()

    def test_returns_absolute_path(self, tmp_path):
        fig, _ = plt.subplots()
        result = _save_and_close(fig, 72, str(tmp_path), "out.png")
        assert os.path.isabs(result)
        assert result.endswith("out.png")


# ===========================================================================
# plot_offtarget_logo  (rendering — just check return type)
# ===========================================================================

class TestPlotOfftargetLogo:
    def test_returns_figure(self, minimal_df):
        matrix = build_offtarget_logo_matrix(minimal_df, guide=None, strip_pam=False)
        fig = plot_offtarget_logo(matrix, title="test", debug=False)
        assert hasattr(fig, "savefig")   # duck-type: it's a Figure
        plt.close("all")

    def test_empty_matrix_raises(self):
        empty = pd.DataFrame(columns=list(LOGO_ALPHABET))
        with pytest.raises(CrispritzReportError):
            plot_offtarget_logo(empty, debug=False)


# ===========================================================================
# plot_annotation_radar  (output-list contract)
# ===========================================================================

class TestPlotAnnotationRadar:
    def test_appends_one_path_per_mm_level(self, minimal_df, tmp_path):
        outputs: list[str] = []
        result = plot_annotation_radar(
            minimal_df,
            annotation_columns=["gene_track"],
            max_mm=1,
            title="test",
            dpi=72,
            outdir=str(tmp_path),
            prefix="rep",
            outputs=outputs,
        )
        # max_mm=1 → levels 0 and 1 → 2 files
        assert len(result) == 2
        assert result is outputs   # same list object (mutated in place)
        plt.close("all")

    def test_output_files_exist(self, minimal_df, tmp_path):
        outputs: list[str] = []
        plot_annotation_radar(
            minimal_df,
            annotation_columns=["gene_track"],
            max_mm=0,
            title="t",
            dpi=72,
            outdir=str(tmp_path),
            prefix="rep",
            outputs=outputs,
        )
        for path in outputs:
            assert os.path.isfile(path)
        plt.close("all")