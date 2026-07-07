"""Partition fan-out and shard merge.

Covers the area "partition fan-out/merge":
  * contig recovery from a partition filename (pure logic),
  * the real C++ k-way merge over synthetic scored shards,
  * the Python fan-out in search_offtargets_tst, with the C++ boundary mocked
    so the orchestration (one task per partition, merge fed every shard) is
    exercised without real .bin files.

Requires the compiled extension (the module imports the C++ core at import).
"""

import pytest

pytest.importorskip("crispritz_plus._ternary_search_tree")

from types import SimpleNamespace

from crispritz_plus.crispritz_cpp import merge_sorted_shards_cpp
from crispritz_plus.search import tst_explorer

_HEADER = "\t".join(
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


def test_contig_from_partition_single_token():
    assert tst_explorer._contig_from_partition("/x/NGG_chr1_1.bin", False) == "chr1"


def test_contig_from_partition_multi_token():
    got = tst_explorer._contig_from_partition("/x/NGG_chr_unplaced_3.bin", False)
    assert got == "chr_unplaced"


def _write_shard(path, rows):
    path.write_text(_HEADER + "\n" + "\n".join(rows) + "\n")


def test_merge_orders_by_edit_distance(tmp_path):
    s1 = tmp_path / "a.shard.tsv"
    s2 = tmp_path / "b.shard.tsv"
    _write_shard(s1, ["chr1\t100\t+\tACGTANGG\tacgtaNGG\t2\tX\t0\t0\t0.10"])
    _write_shard(s2, ["chr1\t200\t+\tACGTANGG\tACGTANGG\t0\tX\t0\t0\t0.99"])
    final = tmp_path / "out.tsv"

    written = merge_sorted_shards_cpp([str(s1), str(s2)], str(final), "edit_distance")
    assert written == 2

    lines = final.read_text().splitlines()
    assert lines[0].split("\t")[0] == "chrom"  # header present
    data = lines[1:]
    # edit_distance sort: the 0-mismatch hit ranks before the 2-mismatch hit
    assert data[0].split("\t")[5] == "0"
    assert data[1].split("\t")[5] == "2"


def test_fan_out_one_task_per_partition(tmp_path, monkeypatch, pam_file, guides_file):
    calls = []
    merged = {}

    def fake_exec(partition, contig, guides, config, pam, upstream, shard, bulge_mode, verbosity):
        calls.append(contig)
        open(shard, "w").close()  # materialise the shard path
        return SimpleNamespace(
            source_path=partition,
            shard_path=shard,
            total_hits=0,
            rows_written=0,
            profiles=[],
        )

    def fake_merge(shard_paths, final_path, sort_mode, *a, **k):
        merged["count"] = len(shard_paths)
        open(final_path, "w").close()
        return 0

    monkeypatch.setattr(tst_explorer, "run_search_executor_cpp", fake_exec)
    monkeypatch.setattr(tst_explorer, "merge_sorted_shards_cpp", fake_merge)
    monkeypatch.setattr(tst_explorer, "write_merged_profiles_cpp", lambda *a, **k: 0)

    outdir = tmp_path / "out"
    outdir.mkdir()
    partitions = [str(tmp_path / "NGG_chr1_1.bin"), str(tmp_path / "NGG_chr2_1.bin")]
    for p in partitions:
        open(p, "w").close()

    tst_explorer.search_offtargets_tst(
        partitions, pam_file, guides_file, 0, 0, 0, str(outdir), 1, 0, False
    )

    assert sorted(calls) == ["chr1", "chr2"]  # fan-out: one task per partition
    assert merged["count"] == 2  # merge fed both shards
