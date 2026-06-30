"""tests/core/test_dna_alphabet.py — unit tests for dna_alphabet.py"""

from __future__ import annotations

import pytest

from crispritz_plus.dna_alphabet import (
    DNA,
    IUPAC,
    RC,
    IUPACTABLE,
    IUPAC_ENCODER,
    reverse_complement,
    dna2rna,
)


# ===========================================================================
# Constants
# ===========================================================================


class TestConstants:
    def test_dna_has_five_symbols(self):
        assert len(DNA) == 5

    def test_dna_contains_n(self):
        assert "N" in DNA

    def test_iupac_extends_dna(self):
        assert all(b in IUPAC for b in DNA)
        assert len(IUPAC) == 15

    def test_rc_covers_all_iupac_uppercase(self):
        for sym in IUPAC:
            assert sym in RC, f"RC missing entry for {sym!r}"

    def test_rc_covers_lowercase(self):
        for sym in IUPAC:
            assert sym.lower() in RC

    def test_iupactable_n_expands_to_acgt(self):
        assert set(IUPACTABLE["N"]) == set("ACGT")

    def test_iupactable_r_expands_to_ag(self):
        assert set(IUPACTABLE["R"]) == {"A", "G"}

    def test_iupac_encoder_ag_and_ga_both_map_to_r(self):
        assert IUPAC_ENCODER["AG"] == "R"
        assert IUPAC_ENCODER["GA"] == "R"

    def test_iupac_encoder_single_bases_are_identity(self):
        for base in "ACGT":
            assert IUPAC_ENCODER[base] == base


# ===========================================================================
# reverse_complement
# ===========================================================================


class TestReverseComplement:
    @pytest.mark.parametrize(
        "seq, expected",
        [
            ("A", "T"),
            ("T", "A"),
            ("C", "G"),
            ("G", "C"),
            ("ACGT", "ACGT"),  # palindrome
            ("AACG", "CGTT"),
            ("TTTT", "AAAA"),
        ],
    )
    def test_basic(self, seq, expected):
        assert reverse_complement(seq) == expected

    def test_iupac_r_complements_to_y(self):
        assert reverse_complement("R") == "Y"

    def test_lowercase_preserved(self):
        assert reverse_complement("acgt") == "acgt"

    def test_double_rc_is_identity(self):
        seq = "ACGTRYMKSW"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_unknown_character_raises(self):
        with pytest.raises(KeyError):
            reverse_complement("Z")


# ===========================================================================
# dna2rna
# ===========================================================================


class TestDna2Rna:
    @pytest.mark.parametrize(
        "dna, rna",
        [
            ("ACGT", "ACGU"),
            ("acgt", "acgu"),
            ("T", "U"),
            ("t", "u"),
            ("ACGN", "ACGN"),  # N unaffected
            ("", ""),
        ],
    )
    def test_conversion(self, dna, rna):
        assert dna2rna(dna) == rna

    def test_iupac_codes_unaffected(self):
        # Only T/t are replaced; other IUPAC ambiguity codes pass through.
        assert dna2rna("RY") == "RY"
