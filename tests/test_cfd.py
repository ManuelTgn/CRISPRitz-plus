"""CFD application: the compute_cfd kernel and the per-shard file scorer.

Covers the area "CFD application". compute_cfd is tested directly with small
synthetic penalty tables (no model pickles needed); the shard scorer is tested
for its file mechanics with compute_cfd and the model loader monkeypatched, so
the test is independent of the CFD tables and of PAM-stripping details.
"""

import pytest

from crispritz_plus.scores.cfd.cfdscore import compute_cfd
from crispritz_plus.scores import shard_scoring


# --- compute_cfd ------------------------------------------------------------

GUIDE = "ACGTACGTACGTACGTACGT"  # 20 nt
PAM = {"GG": 0.9, "AG": 0.2}


def test_perfect_match_is_just_the_pam_factor():
    assert compute_cfd(GUIDE, GUIDE, "GG", {}, PAM) == pytest.approx(0.9)


def test_lowercase_target_is_not_a_mismatch():
    # The target marks mismatches in lowercase; a matched lowercase base must
    # not be scored as a mismatch.
    assert compute_cfd(GUIDE, GUIDE.lower(), "GG", {}, PAM) == pytest.approx(0.9)


def test_single_mismatch_multiplies_penalty():
    target = "ACATACGTACGTACGTACGT"  # position 3: G -> A
    mm = {"rG:dT,3": 0.4}
    assert compute_cfd(GUIDE, target, "GG", mm, PAM) == pytest.approx(0.4 * 0.9)


def test_bulge_gap_is_neutral():
    target = "AC-TACGTACGTACGTACGT"
    assert compute_cfd(GUIDE, target, "GG", {}, PAM) == pytest.approx(0.9)


def test_unknown_pam_dinucleotide_is_neutral():
    assert compute_cfd(GUIDE, GUIDE, "NN", {}, PAM) == pytest.approx(1.0)


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        compute_cfd("ACGT", "ACG", "GG", {}, PAM)


# --- score_shard_file -------------------------------------------------------


@pytest.fixture
def patched_scorer(monkeypatch):
    # Deterministic score, no model pickles; isolates file mechanics.
    monkeypatch.setattr(shard_scoring, "_ensure_models", lambda debug: ({}, {}))
    monkeypatch.setattr(shard_scoring, "compute_cfd", lambda *a, **k: 0.5)


def _write_shard(path, with_header=True):
    # 10-column scored schema; cfd_score starts as the NA sentinel.
    rows = [
        "chr1\t100\t+\tACGTANGG\tacgtaNGG\t1\tX\t0\t0\tNA",
        "chr1\t200\t-\tACGTANGG\tacgtaTGG\t0\tX\t0\t0\tNA",
    ]
    lines = []
    if with_header:
        lines.append(
            "\t".join(
                [
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
            )
        )
    lines += rows
    path.write_text("\n".join(lines) + "\n")


def test_score_shard_fills_score_column(tmp_path, patched_scorer):
    shard = tmp_path / "p.shard.tsv"
    _write_shard(shard)
    scored = shard_scoring.score_shard_file(str(shard))
    assert scored == 2
    out_lines = shard.read_text().splitlines()
    # header preserved
    assert out_lines[0].split("\t")[0] == "chrom"
    # every data row's last column is now the computed score, not NA
    for line in out_lines[1:]:
        assert line.split("\t")[-1] == "0.50"


def test_score_shard_without_header(tmp_path, patched_scorer):
    shard = tmp_path / "p.shard.tsv"
    _write_shard(shard, with_header=False)
    scored = shard_scoring.score_shard_file(str(shard))
    assert scored == 2
    for line in shard.read_text().splitlines():
        assert line.split("\t")[-1] == "0.50"


def test_score_shard_is_in_place(tmp_path, patched_scorer):
    shard = tmp_path / "p.shard.tsv"
    _write_shard(shard)
    shard_scoring.score_shard_file(str(shard))
    # the temp file must not survive the atomic replace
    assert not (tmp_path / "p.shard.tsv.scored.tmp").exists()
