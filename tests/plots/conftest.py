"""
tests/plots/conftest.py
=======================
Bootstrap for the plots test suite.

Architecture note
-----------------
Stubs are injected into ``sys.modules`` **at module-load time** (before
pytest collects test files) so that the top-level imports in
``test_plots.py`` resolve against our controlled stubs rather than the
real package.  A session-scoped autouse fixture restores the original
``sys.modules`` state after all plots tests finish, preventing leakage
into other test sub-trees in the same session.

``logomaker`` is also stubbed here because it is an optional/heavy
install that is not required for unit-testing the pure-logic paths in
``plots.py``.  Rendering tests that do call ``plot_offtarget_logo`` or
``_plot_radar`` receive a ``MagicMock`` logo object and still exercise
the matplotlib figure creation.
"""

from __future__ import annotations

import io
import os
import sys
import types
from unittest.mock import MagicMock

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------


class CrispritzReportError(Exception):
    """Minimal stand-in; raises so pytest.raises() works."""


def _stub_eh(exc_type, msg, exit_code, debug, cause=None):
    raise exc_type(msg) from cause


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_CONF_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_CONF_DIR, "..", ".."))
_PKG_DIR = os.path.join(_PROJECT_ROOT, "src", "crispritz_plus")
_PLOTS_DIR = os.path.join(_PKG_DIR, "plots")


# ---------------------------------------------------------------------------
# Snapshot + inject
# ---------------------------------------------------------------------------

_SNAPSHOT: dict = {
    k: v
    for k, v in sys.modules.items()
    if k.startswith("crispritz_plus") or k == "logomaker"
}


def _install_stubs() -> None:
    def _mod(name: str, **attrs) -> types.ModuleType:
        m = types.ModuleType(name)
        m.__package__ = name
        for k, v in attrs.items():
            setattr(m, k, v)
        return m

    # crispritz_plus (parent – __path__ lets Python find sub-packages on disk)
    pkg = _mod("crispritz_plus")
    pkg.__path__ = [_PKG_DIR]
    pkg.__file__ = os.path.join(_PKG_DIR, "__init__.py")
    sys.modules.setdefault("crispritz_plus", pkg)  # don't overwrite real install

    # crispritz_plus.plots
    plots_pkg = _mod("crispritz_plus.plots")
    plots_pkg.__path__ = [_PLOTS_DIR]
    plots_pkg.__file__ = os.path.join(_PLOTS_DIR, "__init__.py")
    sys.modules.setdefault("crispritz_plus.plots", plots_pkg)

    # exception_handlers
    sys.modules["crispritz_plus.exception_handlers"] = _mod(
        "crispritz_plus.exception_handlers", exception_handler=_stub_eh
    )

    # crispritz_report_errors
    cre = _mod("crispritz_plus.plots.crispritz_report_errors")
    cre.CrispritzReportError = CrispritzReportError
    sys.modules["crispritz_plus.plots.crispritz_report_errors"] = cre

    # logomaker (stub so tests don't need the library installed)
    logo_mock = MagicMock()
    logo_mock.Logo.return_value = MagicMock()
    sys.modules["logomaker"] = logo_mock


_install_stubs()

# Import the module AFTER stubs are in place so relative imports resolve.
import crispritz_plus.plots.plots as _plots  # noqa: E402


# ---------------------------------------------------------------------------
# Session teardown — restore sys.modules
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session", autouse=True)
def _restore_modules():
    yield
    for key in list(sys.modules):
        if key.startswith("crispritz_plus") or key == "logomaker":
            sys.modules.pop(key, None)
    sys.modules.update(_SNAPSHOT)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def exc_cls():
    return CrispritzReportError


@pytest.fixture
def minimal_df() -> pd.DataFrame:
    """Annotated targets DataFrame: two rows, one mismatch, one annotation col."""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "pos": ["100", "200"],
            "strand": ["+", "+"],
            "grna": ["ACGTACGT", "ACGTACGT"],
            "spacer": ["ACGTACGT", "ATGTACGT"],  # row 1 has T mismatch at pos 1
            "mismatches": [0, 1],
            "bulge_type": ["X", "X"],
            "bulge_dna": ["0", "0"],
            "bulge_rna": ["0", "0"],
            "cfd_score": ["1.0", "0.8"],
            "gene_track": ["promoter", "intron"],
        }
    )


@pytest.fixture
def annotated_tsv(tmp_path, minimal_df) -> str:
    """Write minimal_df (mismatches as str) to a TSV and return its path."""
    df = minimal_df.copy()
    df["mismatches"] = df["mismatches"].astype(str)
    p = tmp_path / "targets.annotated.tsv"
    df.to_csv(str(p), sep="\t", index=False)
    return str(p)
