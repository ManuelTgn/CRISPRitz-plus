"""tests/core/test_pam.py — unit tests for pam.py"""

from __future__ import annotations

import pytest

from crispritz_plus.pam import (
    PAM,
    SPCAS9PAM,
    XCAS9PAM,
    CPF1PAM,
    SACAS9PAM,
    CASXPAM,
    SPCAS9,
    XCAS9,
    CPF1,
    SACAS9,
    CASX,
)
from crispritz_plus.crispritz_errors import CrispritzPamError


def _write_pam(tmp_path, line: str) -> str:
    p = tmp_path / "pam.txt"
    p.write_text(line + "\n")
    return str(p)


class TestPamParsing:
    def test_downstream_pam_size(self, pam_file):
        pam = PAM(pam_file, debug=True)
        assert pam.size == 3  # "NGG"

    def test_downstream_guide_size(self, pam_file):
        pam = PAM(pam_file, debug=True)
        assert pam.guide_size == 20  # 23 chars (21N+GG) - 3 PAM

    def test_downstream_not_upstream(self, pam_file):
        assert PAM(pam_file, debug=True).upstream is False

    def test_upstream_pam_flag(self, upstream_pam_file):
        pam = PAM(upstream_pam_file, debug=True)
        assert pam.upstream is True

    def test_pamseq_is_just_pam_portion(self, pam_file):
        pam = PAM(pam_file, debug=True)
        assert pam.pamseq == "NGG"

    def test_upstream_pamseq_extracted(self, upstream_pam_file):
        pam = PAM(upstream_pam_file, debug=True)
        assert pam.pamseq == "TTTN"

    def test_missing_file_raises(self):
        with pytest.raises((CrispritzPamError, Exception)):
            PAM("/no/such/pam.txt", debug=True)


class TestPamCasSystem:
    def test_spcas9_ngg_recognised(self, tmp_path):
        p = _write_pam(tmp_path, "NNNNNNNNNNNNNNNNNNNNGG 3")
        assert PAM(str(p), debug=True).cas_system == SPCAS9

    def test_xcas9_ngk_recognised(self, tmp_path):
        p = _write_pam(tmp_path, "NNNNNNNNNNNNNNNNNNNNNGK 3")
        assert PAM(str(p), debug=True).cas_system == XCAS9

    def test_cpf1_tttn_upstream_recognised(self, tmp_path):
        p = _write_pam(tmp_path, "TTTNNNNNNNNNNNNNNNNNNNNNN -4")
        assert PAM(str(p), debug=True).cas_system == CPF1

    def test_unknown_pam_gives_minus_one(self, tmp_path):
        p = _write_pam(tmp_path, "NNNNNNNNNNNNNNNNNNNNQQQ 3")
        assert PAM(str(p), debug=True).cas_system == -1

    def test_len_returns_pam_size(self, pam_file):
        pam = PAM(pam_file, debug=True)
        assert len(pam) == 3
