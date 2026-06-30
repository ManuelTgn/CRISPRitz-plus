"""tests/core/test_enricher.py — unit tests for enricher.py pure helpers"""

from __future__ import annotations

import os
import pytest

from crispritz_plus.enrichment.enricher import (
    VARIANTGENOMEDIR,
    SNPDIR,
    INDELSDIR,
    _compute_vid,
    _compute_indel_coordinates,
    _skip_variant,
    _extract_af_idx,
    _split_contigs,
    _prepare_output_dir,
)
from crispritz_plus.genome_io import INDELOFFSET
from crispritz_plus.enrichment.enrichment_pair import EnrichPair


# ===========================================================================
# Module constants
# ===========================================================================


class TestConstants:
    def test_variantgenomedir_value(self):
        assert VARIANTGENOMEDIR == "variants_genome"

    def test_snpdir_under_variantgenomedir(self):
        assert VARIANTGENOMEDIR in SNPDIR

    def test_indelsdir_under_variantgenomedir(self):
        assert VARIANTGENOMEDIR in INDELSDIR


# ===========================================================================
# _compute_vid
# ===========================================================================


class TestComputeVid:
    def test_format_contains_all_parts(self):
        vid = _compute_vid("chr1", 101, "A", "T")
        assert "chr1" in vid
        assert "101" in vid
        assert "A" in vid
        assert "T" in vid

    def test_multi_alt(self):
        vid = _compute_vid("chrX", 500, "G", "A,C")
        assert "A,C" in vid


# ===========================================================================
# _compute_indel_coordinates
# ===========================================================================


class TestComputeIndelCoordinates:
    def test_start_is_pos_minus_offset(self):
        start, _ = _compute_indel_coordinates("A", 100)
        assert start == 100 - INDELOFFSET

    def test_stop_accounts_for_ref_length(self):
        _, stop = _compute_indel_coordinates("ATG", 100)
        assert stop == 100 + INDELOFFSET + 3

    def test_single_base_ref(self):
        start, stop = _compute_indel_coordinates("A", 50)
        assert stop - start == 2 * INDELOFFSET + 1


# ===========================================================================
# _skip_variant
# ===========================================================================


class TestSkipVariant:
    def _variant(self, filter_val):
        """Build a minimal 7-column VCF row with the given FILTER."""
        return ["chr1", "100", ".", "A", "T", "60", filter_val]

    def test_pass_variant_not_skipped(self):
        assert _skip_variant("PASS", keep=False) is False

    def test_non_pass_skipped_when_keep_false(self):
        assert _skip_variant("LowQual", keep=False) is True

    def test_non_pass_kept_when_keep_true(self):
        assert _skip_variant("LowQual", keep=True) is False

    def test_dot_filter_skipped_without_keep(self):
        assert _skip_variant(".", keep=False) is True


# ===========================================================================
# _extract_af_idx
# ===========================================================================


class TestExtractAfIdx:
    def test_finds_af_at_correct_position(self):
        info = "DP=30;AF=0.5;MQ=60"
        assert _extract_af_idx(info, debug=False) == 1

    def test_af_first_field(self):
        info = "AF=0.1;DP=20"
        assert _extract_af_idx(info, debug=False) == 0


# ===========================================================================
# _split_contigs
# ===========================================================================


class TestSplitContigs:
    def test_separates_by_vcf_presence(self):
        fmap = {
            "chr1": EnrichPair(False),
            "chr2": EnrichPair(False),
            "chr3": EnrichPair(False),
        }
        fmap["chr1"].fasta = "chr1.fa"
        fmap["chr1"].vcf = "chr1.vcf.gz"
        fmap["chr2"].fasta = "chr2.fa"
        fmap["chr3"].fasta = "chr3.fa"
        fmap["chr3"].vcf = "chr3.vcf.gz"
        with_vcf, without_vcf = _split_contigs(fmap, verbosity=0)
        assert set(with_vcf) == {"chr1", "chr3"}
        assert set(without_vcf) == {"chr2"}

    def test_all_without_vcf(self):
        fmap = {"chrX": EnrichPair(False)}
        fmap["chrX"].fasta = "x.fa"
        with_vcf, without_vcf = _split_contigs(fmap, verbosity=0)
        assert with_vcf == []
        assert "chrX" in without_vcf


# ===========================================================================
# _prepare_output_dir
# ===========================================================================


class TestPrepareOutputDir:
    def test_creates_snp_and_indels_dirs(self, tmp_path):
        snp_dir, indels_dir = _prepare_output_dir(str(tmp_path), verbosity=0)
        assert os.path.isdir(snp_dir)
        assert os.path.isdir(indels_dir)

    def test_returns_absolute_paths(self, tmp_path):
        snp_dir, indels_dir = _prepare_output_dir(str(tmp_path), verbosity=0)
        assert os.path.isabs(snp_dir)
        assert os.path.isabs(indels_dir)
