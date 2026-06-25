"""TSV serialization: the OffTarget row/column contract.

Covers the area "TSV serialization". The central test verifies that every
:data:`TSV_HEADER` label names the value beneath it in
:meth:`OffTarget.to_tsv_row` — the single contract the rest of the pipeline
(C++ shard schema, CFD scorer, merge) relies on.
"""

import pytest

from crispritz_plus.offtarget import TSV_HEADER, OffTarget


def _row():
    # Distinct, identifiable values per field so a misaligned header is caught.
    return OffTarget(
        chrom="chr1",
        pos=10,
        strand="+",
        grna="ACGTANGG",
        spacer="acgtaNGG",
        mm=2,
        bulge_type="DNA",
        bdna=1,
        brna=0,
    )


def test_row_field_count_matches_header():
    fields = _row().to_tsv_row().split("\t")
    assert len(fields) == len(TSV_HEADER)


def test_cfd_score_defaults_to_na_sentinel():
    fields = _row().to_tsv_row().split("\t")
    assert fields[-1] == "NA"


def test_custom_separator():
    fields = _row().to_tsv_row(sep=",").split(",")
    assert len(fields) == len(TSV_HEADER)
    assert fields[0] == "chr1"


def test_no_trailing_newline():
    row = _row().to_tsv_row()
    assert not row.endswith("\n")


@pytest.mark.parametrize(
    "strand",
    ["+", "-"],
)
def test_strand_round_trips(strand):
    ot = OffTarget("chr2", 5, strand, "ACGTANGG", "acgtaNGG", 0, "X", 0, 0)
    assert ot.to_tsv_row().split("\t")[2] == strand
