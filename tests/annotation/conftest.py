# tests/annotation/conftest.py

import pytest
import sys
import types
from unittest.mock import MagicMock
from typing import List


class CrispritzAnnotationError(Exception):
    pass


def _stub_exception_handler(exc_type, msg, exit_code, debug, cause=None):
    raise exc_type(msg) from cause


@pytest.fixture(scope="session", autouse=True)
def _annotation_stubs(tmp_path_factory):
    """
    Inject crispritz_plus stubs for the duration of the annotation test
    session, then restore whatever was in sys.modules before.
    """
    import os

    tmp = tmp_path_factory.mktemp("pkg")
    pkg_dir = tmp / "crispritz_plus"
    ann_dir = pkg_dir / "annotation"
    ann_dir.mkdir(parents=True)

    # Snapshot what exists before we touch anything
    before = {k: v for k, v in sys.modules.items() if k.startswith("crispritz_plus")}

    def _make(name):
        m = types.ModuleType(name)
        m.__package__ = name
        return m

    pkg = _make("crispritz_plus")
    pkg.__path__ = [str(pkg_dir)]
    pkg.__file__ = str(pkg_dir / "__init__.py")
    sys.modules["crispritz_plus"] = pkg

    ann_pkg = _make("crispritz_plus.annotation")
    ann_pkg.__path__ = [str(ann_dir)]
    ann_pkg.__file__ = str(ann_dir / "__init__.py")
    sys.modules["crispritz_plus.annotation"] = ann_pkg

    eh = _make("crispritz_plus.exception_handlers")
    eh.exception_handler = _stub_exception_handler
    sys.modules["crispritz_plus.exception_handlers"] = eh

    ut = _make("crispritz_plus.utils")
    ut.rename_files = MagicMock()
    ut.remove_file = MagicMock()
    sys.modules["crispritz_plus.utils"] = ut

    vb = _make("crispritz_plus.verbosity")
    vb.VERBOSITY_LVL = {0: 0, 1: 1, 2: 2, 3: 3}
    vb.print_verbosity = MagicMock()
    sys.modules["crispritz_plus.verbosity"] = vb

    cae = _make("crispritz_plus.annotation.crispritz_annotation_error")
    cae.CrispritzAnnotationError = CrispritzAnnotationError
    sys.modules["crispritz_plus.annotation.crispritz_annotation_error"] = cae

    yield  # ← all annotation tests run here

    # Restore: remove stubs we added, put back whatever was there before
    for key in list(sys.modules):
        if key.startswith("crispritz_plus"):
            sys.modules.pop(key, None)
    sys.modules.update(before)


@pytest.fixture
def valid_header() -> List[str]:
    """Complete, valid search-output TSV header (all 10 columns)."""
    return [
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
    ]


@pytest.fixture
def bed_file(tmp_path) -> str:
    """Plain three-record BED file intentionally *out of position order*
    (chr2 record appears between the two chr1 records) so that sorting
    tests have something to reorder."""
    p = tmp_path / "test.bed"
    p.write_text(
        "chr1\t300\t400\tgene_B\n"  # chr1 – later position
        "chr2\t50\t150\tgene_C\n"  # chr2 – should sort after chr1
        "chr1\t100\t200\tgene_A\n"  # chr1 – earlier position
    )
    return str(p)


@pytest.fixture
def targets_tsv(tmp_path, valid_header) -> str:
    """Minimal targets TSV: header + one forward-strand data row.

    The data row places a 4-bp spacer (ACGT, no gaps) at 1-based
    position 101 on chr1 forward strand, giving a 0-based interval
    [100, 104).
    """
    header = "\t".join(valid_header)
    row = "chr1\t101\t+\tACGT\tACGT\t0\tX\t0\t0\t1.0"
    p = tmp_path / "results.tsv"
    p.write_text(f"{header}\n{row}\n")
    return str(p)
