"""Shared pytest fixtures and skip markers for the crispritz-plus suite.

The pure-Python units (CFD scoring, TSV serialization) run anywhere. The
orchestration and search-command tests need the compiled extension
``crispritz_plus._ternary_search_tree``; those are gated behind
:data:`requires_cpp` and skip cleanly when the extension has not been built.
"""

import importlib

import pytest


def _cpp_available() -> bool:
    try:
        importlib.import_module("crispritz_plus._ternary_search_tree")
        return True
    except Exception:
        return False


CPP_AVAILABLE = _cpp_available()

#: Skip a test (or module) when the compiled extension is unavailable.
requires_cpp = pytest.mark.skipif(
    not CPP_AVAILABLE,
    reason="compiled extension 'crispritz_plus._ternary_search_tree' not built",
)


# --- Standard SpCas9 fixture: 5 bp guide "ACGTA" + 3 bp "NGG" PAM ------------


@pytest.fixture
def pam_file(tmp_path):
    """A PAM file for a 5 bp guide + NGG PAM (guide-length Ns + PAM, size 3)."""
    path = tmp_path / "pam.txt"
    path.write_text("NNNNNNGG 3\n")
    return str(path)


@pytest.fixture
def guides_file(tmp_path):
    """A guides file: body "ACGTA" plus the 3 'N' PAM placeholders."""
    path = tmp_path / "guides.txt"
    path.write_text("ACGTANNN\n")
    return str(path)
