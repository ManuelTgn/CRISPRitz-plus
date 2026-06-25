"""End-to-end search command.

Covers the area "search command": build a tiny index from a synthetic genome,
run the full search_offtargets_tst pipeline (fan-out -> shard write -> merge),
and assert the targets table is produced with the on-target hit.

Requires the compiled extension. The PAM/guides file formats come from the
shared fixtures (5 bp guide "ACGTA" + NGG), matching the C++ test fixture.
"""

import glob
import os

import pytest

pytest.importorskip("crispritz_plus._ternary_search_tree")

from crispritz_plus.crispritz_cpp import build_tree_cpp
from crispritz_plus.search.tst_explorer import search_offtargets_tst


def test_search_end_to_end_finds_on_target(tmp_path, pam_file, guides_file):
    # Genome embeds the on-target protospacer ACGTA + PAM TGG (matches NGG).
    genome = "TTTT" + "ACGTATGG" + "TTTT"
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    # build_tree_cpp(sequence, contig, pam, pam_length, pam_size, upstream,
    #                outdir, max_bulges, threads)
    build_tree_cpp(genome, "chr1", "NGG", 8, 3, False, str(index_dir), 0, 1)

    partitions = glob.glob(os.path.join(str(index_dir), "*.bin"))
    assert partitions, "index build produced no .bin partition"

    outdir = tmp_path / "out"
    outdir.mkdir()
    search_offtargets_tst(
        partitions, pam_file, guides_file, 0, 0, 0, str(outdir), 1, 0, False
    )

    final = outdir / "guides.targets.tsv"
    assert final.exists(), "search produced no targets table"

    lines = final.read_text().splitlines()
    data = [ln for ln in lines if ln and not ln.startswith("chrom")]
    assert data, "expected at least the on-target hit"
    assert any(ln.split("\t")[0] == "chr1" for ln in data)
