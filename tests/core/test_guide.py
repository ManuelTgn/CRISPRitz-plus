"""tests/core/test_guide.py — unit tests for guide.py"""

from __future__ import annotations
import pytest
from crispritz_plus.guide import Guide, GuideList, _validate_guide_sequence
from crispritz_plus.pam import PAM
from crispritz_plus.crispritz_errors import CrispritzGuideError

# Guide format rules (enforced by _validate_guide_sequence):
#   - Downstream PAM (e.g. NGG): guide ends with N×pam_size  →  "BODY" + "NNN"
#   - Upstream   PAM (e.g. TTTN): guide starts with N×pam_size → "NNNN" + "BODY"


class TestGuide:
    def test_sequence_strips_downstream_pam(self, pam_file):
        pam = PAM(pam_file, debug=True)
        g = Guide("ACGTACGTACGTACGTACGTNNN", pam, debug=True)
        assert g.sequence == "ACGTACGTACGTACGT ACGT".replace(" ", "")

    def test_len_is_guide_body_length(self, pam_file):
        pam = PAM(pam_file, debug=True)
        g = Guide("ACGTACGTACGTACGTACGTNNN", pam, debug=True)
        assert len(g) == 20

    def test_reverse_is_reverse_complement(self, pam_file):
        pam = PAM(pam_file, debug=True)
        g = Guide("ACGTACGTACGTACGTACGTNNN", pam, debug=True)
        from crispritz_plus.dna_alphabet import reverse_complement

        assert g.reverse == reverse_complement(g.sequence)

    def test_sequence_strips_upstream_pam(self, upstream_pam_file):
        pam = PAM(upstream_pam_file, debug=True)
        # 4 Ns (PAM placeholder) + 20-nt body
        body = "ACGTACGTACGTACGTACGT"
        g = Guide("NNNN" + body, pam, debug=True)
        assert g.sequence == body

    def test_bad_pam_placeholder_raises(self, pam_file):
        pam = PAM(pam_file, debug=True)
        # Last 3 chars are "GGG" instead of "NNN"
        with pytest.raises(CrispritzGuideError):
            Guide("ACGTACGTACGTACGTACGTGGG", pam, debug=True)


class TestValidateGuideSequence:
    def test_valid_downstream_returns_body(self):
        seq = _validate_guide_sequence("ACGTACGTACGTACGTACGTNNN", 3, False, True)
        assert seq == "ACGTACGTACGTACGTACGT"

    def test_valid_upstream_returns_body(self):
        seq = _validate_guide_sequence("NNNACGTACGTACGTACGTACGT", 3, True, True)
        assert seq == "ACGTACGTACGTACGTACGT"

    def test_non_n_in_downstream_pam_slot_raises(self):
        with pytest.raises(CrispritzGuideError):
            _validate_guide_sequence("ACGTACGTACGTACGTACGTGGG", 3, False, True)

    def test_non_n_in_upstream_pam_slot_raises(self):
        with pytest.raises(CrispritzGuideError):
            _validate_guide_sequence("TTTACGTACGTACGTACGTACGT", 3, True, True)


class TestGuideList:
    def test_reads_two_guides(self, guide_file, pam_file):
        pam = PAM(pam_file, debug=True)
        gl = GuideList(guide_file, pam, debug=True)
        assert len(gl.guides) == 2

    def test_guides_property_returns_list(self, guide_file, pam_file):
        pam = PAM(pam_file, debug=True)
        assert isinstance(GuideList(guide_file, pam, debug=True).guides, list)

    def test_empty_guide_file_raises(self, tmp_path, pam_file):
        pam = PAM(pam_file, debug=True)
        p = tmp_path / "empty.txt"
        p.write_text("")
        with pytest.raises(CrispritzGuideError):
            GuideList(str(p), pam, debug=True)

    def test_mismatched_lengths_raises(self, tmp_path, pam_file):
        pam = PAM(pam_file, debug=True)
        p = tmp_path / "bad.txt"
        # Second guide has fewer body bases → different total length from first
        p.write_text("ACGTACGTACGTACGTACGTNNN\nACGTNNN\n")
        with pytest.raises(CrispritzGuideError):
            GuideList(str(p), pam, debug=True)

    def test_missing_file_raises(self, pam_file):
        pam = PAM(pam_file, debug=True)
        with pytest.raises((CrispritzGuideError, Exception)):
            GuideList("/no/such/guides.txt", pam, debug=True)
