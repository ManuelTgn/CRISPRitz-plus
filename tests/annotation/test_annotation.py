"""
tests/annotation/test_annotation.py
=====================================
Unit tests for ``crispritz_plus.annotation.annotation``.

Coverage map
------------
* Module constants          - :class:`TestModuleConstants`
* ``BedRecord`` dataclass   - :class:`TestBedRecord`
* ``OffTarget`` dataclass   - :class:`TestOffTarget`
* ``_create_targets_ann``   - :class:`TestCreateTargetsAnn`
* ``_validate_targets_header`` - :class:`TestValidateTargetsHeader`
* ``_resolve_annotation_names`` - :class:`TestResolveAnnotationNames`
* ``_is_bgzipped``          - :class:`TestIsBgzipped`
* ``_read_bed``             - :class:`TestReadBed`
* ``_sort_bed``             - :class:`TestSortBed`
* ``_genomic_footprint``    - :class:`TestGenomicFootprint`
* ``_compute_target_footprint`` - :class:`TestComputeTargetFootprint`
* ``_overlapping_features`` - :class:`TestOverlappingFeatures`
* ``_update_header``        - :class:`TestUpdateHeader`
* ``_retrieve_offtarget_position`` - :class:`TestRetrieveOfftargetPosition`
* ``_read_offtarget``       - :class:`TestReadOfftarget`
* ``_prepare_bed_track``    - :class:`TestPrepareBedTrack`
* ``_annotate_stream``      - :class:`TestAnnotateStream`

Dependency model
----------------
``conftest.py`` (in this directory) replaces all internal crispritz_plus
imports with stubs before this module is collected.  ``pysam`` is used
as the real library; individual tests mock specific pysam symbols where
real bgzip/tabix binaries are not required.
"""

from __future__ import annotations

import io
import os
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest

# Module under test - importable after conftest has patched sys.modules.
import crispritz_plus.annotation.annotation as ann
from crispritz_plus.annotation.annotation import (
    # Constants
    SEARCH_OUTPUT_HEADER,
    _BED_SKIP_PREFIXES,
    _CORE_COLUMNS,
    _NO_OVERLAP,
    _OUTPUT_SUFFIX,
    # Data classes
    BedRecord,
    OffTarget,
    # Helpers - output path
    _create_targets_ann,
    # Helpers - header
    _validate_targets_header,
    _resolve_annotation_names,
    _update_header,
    # Helpers - BED preparation
    _is_bgzipped,
    _read_bed,
    _sort_bed,
    # Helpers - coordinate arithmetic
    _genomic_footprint,
    _compute_target_footprint,
    # Helpers - per-site annotation
    _overlapping_features,
    _retrieve_offtarget_position,
    _read_offtarget,
    # Pipeline entry points
    _prepare_bed_track,
    _annotate_stream,
)
from crispritz_plus.annotation.crispritz_annotation_error import (
    CrispritzAnnotationError,
)


# ===========================================================================
# Helpers shared across test classes
# ===========================================================================

#: Full column-index mapping matching the 10-column valid header fixture.
_FULL_COL_IDX: dict[str, int] = {
    "chrom": 0,
    "pos": 1,
    "strand": 2,
    "grna": 3,
    "spacer": 4,
    "mismatches": 5,
    "bulge_type": 6,
    "bulge_dna": 7,
    "bulge_rna": 8,
    "cfd_score": 9,
}


def _mock_tabix(rows: list[str]) -> MagicMock:
    """Return a ``pysam.TabixFile`` mock whose ``fetch`` yields *rows*."""
    t = MagicMock()
    t.fetch.return_value = iter(rows)
    return t


def _make_fields(
    chrom: str = "chr1",
    pos: str = "101",
    strand: str = "+",
    grna: str = "ACGT",
    spacer: str = "ACGT",
) -> list[str]:
    """Return a 10-field targets-TSV data row as a list of strings."""
    return [chrom, pos, strand, grna, spacer, "0", "X", "0", "0", "1.0"]


# ===========================================================================
# Module-level constants
# ===========================================================================


class TestModuleConstants:
    def test_search_output_header_is_list(self):
        assert isinstance(SEARCH_OUTPUT_HEADER, list)

    def test_search_output_header_has_ten_columns(self):
        assert len(SEARCH_OUTPUT_HEADER) == 10

    def test_first_column_is_chrom(self):
        assert SEARCH_OUTPUT_HEADER[0] == "chrom"

    def test_core_columns_first_nine_match_header(self):
        assert SEARCH_OUTPUT_HEADER[:_CORE_COLUMNS] == [
            "chrom",
            "pos",
            "strand",
            "grna",
            "spacer",
            "mismatches",
            "bulge_type",
            "bulge_dna",
            "bulge_rna",
        ]

    def test_core_columns_is_nine(self):
        assert _CORE_COLUMNS == 9

    def test_no_overlap_sentinel_is_na(self):
        assert _NO_OVERLAP == "NA"

    def test_output_suffix_value(self):
        assert _OUTPUT_SUFFIX == ".annotated.tsv"

    def test_bed_skip_prefixes_comment(self):
        assert "#" in _BED_SKIP_PREFIXES

    def test_bed_skip_prefixes_track(self):
        assert "track" in _BED_SKIP_PREFIXES

    def test_bed_skip_prefixes_browser(self):
        assert "browser" in _BED_SKIP_PREFIXES


# ===========================================================================
# BedRecord
# ===========================================================================


class TestBedRecord:
    def test_instantiation_and_field_access(self):
        rec = BedRecord(chrom="chr1", pos=100, line="chr1\t100\t200\tgene_A")
        assert rec.chrom == "chr1"
        assert rec.pos == 100
        assert rec.line == "chr1\t100\t200\tgene_A"

    def test_dataclass_equality(self):
        r1 = BedRecord("chr1", 100, "chr1\t100\t200\tgene_A")
        r2 = BedRecord("chr1", 100, "chr1\t100\t200\tgene_A")
        assert r1 == r2

    def test_inequality_when_pos_differs(self):
        r1 = BedRecord("chr1", 100, "chr1\t100\t200\tgene_A")
        r2 = BedRecord("chr1", 200, "chr1\t200\t300\tgene_A")
        assert r1 != r2


# ===========================================================================
# OffTarget
# ===========================================================================


class TestOffTarget:
    def test_forward_strand_instantiation(self):
        ot = OffTarget(chrom="chr1", start=100, end=120, strand="+")
        assert ot.chrom == "chr1"
        assert ot.start == 100
        assert ot.end == 120
        assert ot.strand == "+"

    def test_reverse_strand_instantiation(self):
        ot = OffTarget(chrom="chrX", start=50, end=70, strand="-")
        assert ot.strand == "-"

    def test_dataclass_equality(self):
        a = OffTarget("chr2", 300, 320, "+")
        b = OffTarget("chr2", 300, 320, "+")
        assert a == b

    def test_inequality_when_strand_differs(self):
        a = OffTarget("chr1", 100, 120, "+")
        b = OffTarget("chr1", 100, 120, "-")
        assert a != b


# ===========================================================================
# _create_targets_ann
# ===========================================================================


class TestCreateTargetsAnn:
    def test_basic_stem_and_suffix(self, tmp_path: Path):
        result = _create_targets_ann("results.tsv", str(tmp_path))
        assert result == os.path.join(
            os.path.abspath(str(tmp_path)), "results.annotated.tsv"
        )

    def test_only_last_extension_stripped(self, tmp_path: Path):
        result = _create_targets_ann("targets.search.tsv", str(tmp_path))
        assert os.path.basename(result) == "targets.search.annotated.tsv"

    def test_outdir_resolved_to_absolute(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        result = _create_targets_ann("results.tsv", ".")
        assert os.path.isabs(result)

    def test_basename_used_not_full_input_path(self, tmp_path: Path):
        result = _create_targets_ann("/deep/path/to/my_results.tsv", str(tmp_path))
        assert os.path.basename(result) == "my_results.annotated.tsv"


# ===========================================================================
# _validate_targets_header
# ===========================================================================


class TestValidateTargetsHeader:
    def test_valid_header_returns_full_index_mapping(self, valid_header: list[str]):
        col_idx = _validate_targets_header(valid_header, debug=False)
        assert col_idx["chrom"] == 0
        assert col_idx["pos"] == 1
        assert col_idx["strand"] == 2
        assert col_idx["cfd_score"] == 9

    def test_extra_trailing_columns_are_accepted(self, valid_header: list[str]):
        extended = valid_header + ["extra_col"]
        col_idx = _validate_targets_header(extended, debug=False)
        assert col_idx["extra_col"] == 10

    def test_mapping_covers_all_supplied_columns(self, valid_header: list[str]):
        col_idx = _validate_targets_header(valid_header, debug=False)
        assert set(col_idx.keys()) == set(valid_header)


# ===========================================================================
# _resolve_annotation_names
# ===========================================================================


class TestResolveAnnotationNames:
    def test_none_generates_annotation_n_names(self):
        names = _resolve_annotation_names(
            ["a.bed", "b.bed", "c.bed"], None, debug=False
        )
        assert names == ["annotation1", "annotation2", "annotation3"]

    def test_none_single_track_produces_annotation1(self):
        names = _resolve_annotation_names(["a.bed"], None, debug=False)
        assert names == ["annotation1"]

    def test_none_empty_list_returns_empty(self):
        names = _resolve_annotation_names([], None, debug=False)
        assert names == []

    def test_explicit_names_returned_unchanged(self):
        explicit = ["genes", "repeats"]
        names = _resolve_annotation_names(["a.bed", "b.bed"], explicit, debug=False)
        assert names == explicit


# ===========================================================================
# _is_bgzipped
# ===========================================================================


class TestIsBgzipped:
    @pytest.mark.parametrize(
        "path, expected",
        [
            ("annotations.bed.gz", True),
            ("file.gz", True),
            ("file.bed", False),
            ("file.tsv", False),
            ("", False),
            ("file.GZ", False),  # case-sensitive: only lowercase .gz
        ],
    )
    def test_extension_detection(self, path: str, expected: bool):
        assert _is_bgzipped(path) is expected


# ===========================================================================
# _read_bed
# ===========================================================================


class TestReadBed:
    def test_parses_all_valid_records(self, bed_file: str):
        records = _read_bed(bed_file, debug=False)
        assert len(records) == 3

    def test_correct_chrom_and_pos_parsed(self, bed_file: str):
        records = _read_bed(bed_file, debug=False)
        chroms = {r.chrom for r in records}
        positions = {r.pos for r in records}
        assert chroms == {"chr1", "chr2"}
        assert positions == {300, 50, 100}

    def test_line_preserved_verbatim(self, tmp_path: Path):
        raw = "chr1\t100\t200\tgene_A"
        p = tmp_path / "single.bed"
        p.write_text(raw + "\n")
        records = _read_bed(str(p), debug=False)
        assert len(records) == 1
        assert records[0].line == raw

    def test_skips_hash_comment_lines(self, tmp_path: Path):
        p = tmp_path / "f.bed"
        p.write_text("# comment\nchr1\t10\t20\tgene\n")
        assert len(_read_bed(str(p), debug=False)) == 1

    def test_skips_track_directive(self, tmp_path: Path):
        p = tmp_path / "f.bed"
        p.write_text("track name=foo\nchr1\t10\t20\tgene\n")
        assert len(_read_bed(str(p), debug=False)) == 1

    def test_skips_browser_directive(self, tmp_path: Path):
        p = tmp_path / "f.bed"
        p.write_text("browser position chr1:1-1000\nchr1\t10\t20\tgene\n")
        assert len(_read_bed(str(p), debug=False)) == 1

    def test_skips_blank_lines(self, tmp_path: Path):
        p = tmp_path / "f.bed"
        p.write_text("\nchr1\t10\t20\tgene\n\n")
        assert len(_read_bed(str(p), debug=False)) == 1

    def test_skips_records_with_fewer_than_four_fields(self, tmp_path: Path):
        p = tmp_path / "f.bed"
        p.write_text(
            "chr1\t10\t20\n"  # 3 fields - skipped
            "chr1\t30\t40\tgene\n"  # 4 fields - kept
        )
        assert len(_read_bed(str(p), debug=False)) == 1

    def test_returns_records_in_file_order(self, tmp_path: Path):
        p = tmp_path / "f.bed"
        p.write_text("chrZ\t900\t1000\tZ_gene\n" "chrA\t1\t100\tA_gene\n")
        records = _read_bed(str(p), debug=False)
        assert [r.chrom for r in records] == ["chrZ", "chrA"]


# ===========================================================================
# _sort_bed
# ===========================================================================


class TestSortBed:
    def test_returns_sorted_file_path(self, bed_file: str):
        sorted_path = _sort_bed(bed_file, debug=False)
        assert sorted_path == f"{bed_file}.sorted.bed"

    def test_sorted_file_is_created_on_disk(self, bed_file: str):
        sorted_path = _sort_bed(bed_file, debug=False)
        assert os.path.isfile(sorted_path)

    def test_chromosomes_in_lexicographic_order(self, bed_file: str):
        sorted_path = _sort_bed(bed_file, debug=False)
        with open(sorted_path) as fh:
            chroms = [ln.split("\t")[0] for ln in fh if ln.strip()]
        assert chroms == sorted(chroms)

    def test_positions_sorted_within_chromosome(self, bed_file: str):
        sorted_path = _sort_bed(bed_file, debug=False)
        with open(sorted_path) as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        chr1_pos = [int(ln.split("\t")[1]) for ln in lines if ln.startswith("chr1")]
        assert chr1_pos == sorted(chr1_pos)

    def test_all_records_present_after_sort(self, bed_file: str):
        sorted_path = _sort_bed(bed_file, debug=False)
        with open(sorted_path) as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        assert len(lines) == 3


# ===========================================================================
# _genomic_footprint
# ===========================================================================


class TestGenomicFootprint:
    @pytest.mark.parametrize(
        "spacer, expected",
        [
            ("ACGT", 4),  # no gaps
            ("AC-GT", 4),  # one internal gap
            ("A--T", 2),  # two gaps
            ("----", 0),  # all gaps
            ("", 0),  # empty spacer
            ("NNNNNN", 6),  # non-standard bases, no gaps
            ("-A-C-", 2),  # gaps interspersed
        ],
    )
    def test_footprint(self, spacer: str, expected: int):
        assert _genomic_footprint(spacer) == expected


# ===========================================================================
# _compute_target_footprint
# ===========================================================================


class TestComputeTargetFootprint:
    def test_forward_strand_basic(self):
        start, end = _compute_target_footprint(101, "ACGT", "+")
        # 1-based 101 → 0-based 100; footprint 4 → [100, 104)
        assert start == 100
        assert end == 104

    def test_reverse_strand_basic(self):
        start, end = _compute_target_footprint(101, "ACGT", "-")
        # pos0=100; interval [100-4, 100) = [96, 100)
        assert start == 96
        assert end == 100

    def test_forward_with_gapped_spacer(self):
        # "AC-GT" has footprint 4 - same interval as "ACGT"
        start, end = _compute_target_footprint(101, "AC-GT", "+")
        assert end - start == 4

    def test_reverse_with_gapped_spacer(self):
        start, end = _compute_target_footprint(101, "AC-GT", "-")
        assert end - start == 4

    def test_zero_footprint_forward(self):
        start, end = _compute_target_footprint(101, "", "+")
        assert start == end == 100

    def test_zero_footprint_reverse(self):
        start, end = _compute_target_footprint(101, "", "-")
        assert start == end == 100

    def test_interval_is_half_open(self):
        # [start, end) must have length == footprint
        spacer = "ACGTACGT"  # 8 bp
        start, end = _compute_target_footprint(201, spacer, "+")
        assert end - start == 8


# ===========================================================================
# _overlapping_features
# ===========================================================================


class TestOverlappingFeatures:
    def test_single_overlapping_feature(self):
        tabix = _mock_tabix(["chr1\t90\t110\tpromoter"])
        ot = OffTarget("chr1", 100, 104, "+")
        assert _overlapping_features(tabix, ot) == ["promoter"]

    def test_no_overlap_returns_empty_list(self):
        tabix = _mock_tabix([])
        ot = OffTarget("chr1", 100, 104, "+")
        assert _overlapping_features(tabix, ot) == []

    def test_missing_contig_value_error_returns_empty(self):
        tabix = MagicMock()
        tabix.fetch.side_effect = ValueError("chromosome not in index")
        ot = OffTarget("chrUNKNOWN", 0, 10, "+")
        assert _overlapping_features(tabix, ot) == []

    def test_duplicate_labels_deduplicated(self):
        rows = [
            "chr1\t90\t110\tgene_A",
            "chr1\t95\t115\tgene_A",  # same label
        ]
        labels = _overlapping_features(
            _mock_tabix(rows), OffTarget("chr1", 100, 104, "+")
        )
        assert labels == ["gene_A"]

    def test_multiple_distinct_labels_all_returned(self):
        rows = [
            "chr1\t90\t110\tgene_A",
            "chr1\t95\t115\tgene_B",
        ]
        labels = _overlapping_features(
            _mock_tabix(rows), OffTarget("chr1", 100, 104, "+")
        )
        assert set(labels) == {"gene_A", "gene_B"}
        assert len(labels) == 2

    def test_first_occurrence_order_preserved(self):
        # gene_B appears first; gene_A second; second gene_B is a duplicate
        rows = [
            "chr1\t90\t110\tgene_B",
            "chr1\t92\t112\tgene_A",
            "chr1\t94\t114\tgene_B",  # duplicate
        ]
        labels = _overlapping_features(
            _mock_tabix(rows), OffTarget("chr1", 100, 104, "+")
        )
        assert labels == ["gene_B", "gene_A"]

    def test_fetch_called_with_correct_coordinates(self):
        tabix = _mock_tabix([])
        ot = OffTarget("chrX", 500, 520, "-")
        _overlapping_features(tabix, ot)
        tabix.fetch.assert_called_once_with("chrX", 500, 520)


# ===========================================================================
# _update_header
# ===========================================================================


class TestUpdateHeader:
    def test_new_column_names_appended_to_written_header(self, valid_header: list[str]):
        header_line = "\t".join(valid_header)
        fout = io.StringIO()
        _update_header(
            header_line, fout, ["track1", "track2"], "dummy.tsv", debug=False
        )
        written = fout.getvalue().rstrip("\n")
        assert written.endswith("track1\ttrack2")

    def test_original_columns_preserved_in_written_header(
        self, valid_header: list[str]
    ):
        header_line = "\t".join(valid_header)
        fout = io.StringIO()
        _update_header(header_line, fout, ["extra"], "dummy.tsv", debug=False)
        written = fout.getvalue()
        for col in valid_header:
            assert col in written

    def test_returns_correct_column_index_mapping(self, valid_header: list[str]):
        header_line = "\t".join(valid_header)
        fout = io.StringIO()
        col_idx = _update_header(header_line, fout, ["x"], "dummy.tsv", debug=False)
        assert col_idx["chrom"] == 0
        assert col_idx["pos"] == 1


# ===========================================================================
# _retrieve_offtarget_position
# ===========================================================================


class TestRetrieveOfftargetPosition:
    def test_valid_integer_string_parsed(self):
        fields = ["chr1", "101", "+", "ACGT", "ACGT"]
        assert (
            _retrieve_offtarget_position(fields, idx=1, line_no=2, debug=False) == 101
        )

    def test_position_one_is_valid(self):
        fields = ["chr1", "1"]
        assert _retrieve_offtarget_position(fields, idx=1, line_no=1, debug=False) == 1


# ===========================================================================
# _read_offtarget
# ===========================================================================


class TestReadOfftarget:
    def test_forward_strand_interval(self):
        # pos1=101, spacer="ACGT" (4 bp), strand=+ → [100, 104)
        ot = _read_offtarget(
            _make_fields(pos="101", strand="+", spacer="ACGT"),
            _FULL_COL_IDX,
            line_no=2,
            debug=False,
        )
        assert ot.start == 100
        assert ot.end == 104

    def test_reverse_strand_interval(self):
        # pos1=101, spacer="ACGT" (4 bp), strand=- → [96, 100)
        ot = _read_offtarget(
            _make_fields(pos="101", strand="-", spacer="ACGT"),
            _FULL_COL_IDX,
            line_no=2,
            debug=False,
        )
        assert ot.start == 96
        assert ot.end == 100

    def test_chrom_extracted_correctly(self):
        ot = _read_offtarget(
            _make_fields(chrom="chr7"), _FULL_COL_IDX, line_no=2, debug=False
        )
        assert ot.chrom == "chr7"

    def test_strand_stored_on_offtarget(self):
        ot = _read_offtarget(
            _make_fields(strand="-"), _FULL_COL_IDX, line_no=2, debug=False
        )
        assert ot.strand == "-"

    def test_gapped_spacer_reduces_footprint(self):
        # "AC-GT" has footprint 4 (same as "ACGT") so interval width stays 4
        ot = _read_offtarget(
            _make_fields(pos="101", strand="+", spacer="AC-GT"),
            _FULL_COL_IDX,
            line_no=2,
            debug=False,
        )
        assert ot.end - ot.start == 4

    def test_returns_offtarget_instance(self):
        ot = _read_offtarget(_make_fields(), _FULL_COL_IDX, line_no=2, debug=False)
        assert isinstance(ot, OffTarget)


# ===========================================================================
# _prepare_bed_track
# ===========================================================================


class TestPrepareBedTrack:
    def test_already_bgzipped_and_indexed_passthrough(self, tmp_path: Path):
        gz = tmp_path / "track.bed.gz"
        tbi = tmp_path / "track.bed.gz.tbi"
        gz.write_bytes(b"")
        tbi.write_bytes(b"")
        result = _prepare_bed_track(str(gz), debug=False)
        assert result == str(gz)


# ===========================================================================
# _annotate_stream
# ===========================================================================


class TestAnnotateStream:
    def test_augmented_header_written_to_output(self, targets_tsv: str):
        fout = io.StringIO()
        _annotate_stream(
            targets_tsv,
            fout,
            tabix_handles=[_mock_tabix([])],
            column_names=["my_track"],
            verbosity=0,
            debug=False,
        )
        header = fout.getvalue().split("\n")[0]
        assert "my_track" in header

    def test_all_original_header_columns_preserved(
        self, targets_tsv: str, valid_header: list[str]
    ):
        fout = io.StringIO()
        _annotate_stream(
            targets_tsv, fout, [_mock_tabix([])], ["t"], verbosity=0, debug=False
        )
        header = fout.getvalue().split("\n")[0]
        for col in valid_header:
            assert col in header

    def test_overlapping_feature_written_to_data_row(self, targets_tsv: str):
        tabix = _mock_tabix(["chr1\t90\t110\tpromoter"])
        fout = io.StringIO()
        _annotate_stream(
            targets_tsv, fout, [tabix], ["genes"], verbosity=0, debug=False
        )
        data_line = fout.getvalue().strip().split("\n")[1]
        assert "promoter" in data_line

    def test_no_overlap_writes_na_sentinel(self, targets_tsv: str):
        fout = io.StringIO()
        _annotate_stream(
            targets_tsv, fout, [_mock_tabix([])], ["genes"], verbosity=0, debug=False
        )
        data_line = fout.getvalue().strip().split("\n")[1]
        assert data_line.endswith(_NO_OVERLAP)

    def test_two_tracks_produce_two_annotation_cells(self, targets_tsv: str):
        t1 = _mock_tabix(["chr1\t90\t110\texon"])
        t2 = _mock_tabix([])
        fout = io.StringIO()
        _annotate_stream(
            targets_tsv, fout, [t1, t2], ["t1", "t2"], verbosity=0, debug=False
        )
        data_cells = fout.getvalue().strip().split("\n")[1].split("\t")
        assert data_cells[-2] == "exon"
        assert data_cells[-1] == _NO_OVERLAP

    def test_returns_correct_annotated_count(self, targets_tsv: str):
        fout = io.StringIO()
        count = _annotate_stream(
            targets_tsv, fout, [_mock_tabix([])], ["t"], verbosity=0, debug=False
        )
        assert count == 1

    def test_blank_data_lines_not_counted(
        self, tmp_path: Path, valid_header: list[str]
    ):
        header = "\t".join(valid_header)
        row = "chr1\t101\t+\tACGT\tACGT\t0\tX\t0\t0\t1.0"
        p = tmp_path / "blanks.tsv"
        p.write_text(f"{header}\n{row}\n\n\n")  # extra blank lines
        fout = io.StringIO()
        count = _annotate_stream(
            str(p), fout, [_mock_tabix([])], ["t"], verbosity=0, debug=False
        )
        assert count == 1

    def test_comma_joined_multiple_features_in_cell(self, targets_tsv: str):
        # Two distinct features overlap → joined by comma
        rows = [
            "chr1\t90\t110\tgene_A",
            "chr1\t95\t115\tgene_B",
        ]
        fout = io.StringIO()
        _annotate_stream(
            targets_tsv, fout, [_mock_tabix(rows)], ["track"], verbosity=0, debug=False
        )
        data_cells = fout.getvalue().strip().split("\n")[1].split("\t")
        last_cell = data_cells[-1]
        assert "gene_A" in last_cell
        assert "gene_B" in last_cell
        assert "," in last_cell

    def test_tabix_fetch_called_with_derived_coordinates(self, targets_tsv: str):
        # pos=101, spacer=ACGT, strand=+ → fetch("chr1", 100, 104)
        tabix = _mock_tabix([])
        fout = io.StringIO()
        _annotate_stream(targets_tsv, fout, [tabix], ["t"], verbosity=0, debug=False)
        tabix.fetch.assert_called_once_with("chr1", 100, 104)
