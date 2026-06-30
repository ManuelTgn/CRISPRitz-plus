"""tests/core/conftest.py — shared fixtures for all crispritz_plus unit tests."""

from __future__ import annotations

import os
import sys
import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)


# -------------------------------------------------------------------------
# PAM files
# -------------------------------------------------------------------------


@pytest.fixture
def pam_file(tmp_path):
    """SpCas9 NGG PAM: 20-nt guide body + 3-nt PAM. guide_size = 20."""
    # 21 Ns + GG  = 23 chars; last 3 = NGG → pamseq="NGG", guide_size=20
    p = tmp_path / "pam.txt"
    p.write_text("NNNNNNNNNNNNNNNNNNNNNGG 3\n")
    return str(p)


@pytest.fixture
def upstream_pam_file(tmp_path):
    """Cpf1 TTTN PAM (upstream): 4-nt PAM + 20-nt guide body. guide_size = 20."""
    # TTTN + 20 Ns = 24 chars; first 4 = TTTN → pamseq="TTTN", guide_size=20
    p = tmp_path / "pam_up.txt"
    p.write_text("TTTNNNNNNNNNNNNNNNNNNNNNN -4\n")
    return str(p)


# -------------------------------------------------------------------------
# FASTA / guide files
# -------------------------------------------------------------------------


@pytest.fixture
def fasta_file(tmp_path):
    """Single-record FASTA, 100 bp."""
    seq = "ACGT" * 25
    p = tmp_path / "chr1.fa"
    p.write_text(f">chr1\n{seq}\n")
    return str(p)


@pytest.fixture
def guide_file(tmp_path):
    """Two SpCas9 guides: 20-nt body + NNN PAM placeholder (NOT the actual NGG)."""
    # PAM placeholder = NNN (all N), because _validate_guide_sequence demands it
    p = tmp_path / "guides.txt"
    p.write_text("ACGTACGTACGTACGTACGTNNN\nTTTTAAAACCCCGGGGAAAANNN\n")
    return str(p)
