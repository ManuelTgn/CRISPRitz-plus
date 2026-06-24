"""Config construction: SearchConfiguration and the enum string parsers.

Covers the area "config construction". Requires the compiled extension because
SearchConfiguration and the *_from_string parsers live in the C++ core.
"""

import pytest

pytest.importorskip("crispritz_plus._ternary_search_tree")

from crispritz_plus.crispritz_cpp import make_search_configuration
from crispritz_plus.crispritz_cpp.bulge_mode import BulgeMode
from crispritz_plus.crispritz_cpp.output_mode import OutputMode
from crispritz_plus.crispritz_cpp.sort_mode import SortMode


def test_make_configuration_stores_budgets():
    cfg = make_search_configuration(3, 1, 1, 4, output_mode="both").native
    assert cfg.max_mismatches == 3
    assert cfg.max_bulges_dna == 1
    assert cfg.max_bulges_rna == 1
    assert cfg.threads == 4


def test_max_total_edits_is_sum_of_budgets():
    cfg = make_search_configuration(3, 2, 1, 1, output_mode="both").native
    assert cfg.max_total_edits == 3 + 2 + 1


@pytest.mark.parametrize(
    "mode,write_targets,write_profile",
    [
        ("both", True, True),
        ("targets", True, False),
        ("profile", False, True),
    ],
)
def test_output_mode_predicates(mode, write_targets, write_profile):
    cfg = make_search_configuration(0, 0, 0, 1, output_mode=mode).native
    assert cfg.write_targets == write_targets
    assert cfg.write_profile == write_profile


def test_bulge_mode_from_string_round_trip():
    assert BulgeMode.from_string("mixed") == BulgeMode.MixedBulges
    assert BulgeMode.from_string("single") == BulgeMode.SingleBulgeType


def test_bulge_mode_invalid_token_raises():
    with pytest.raises(Exception):
        BulgeMode.from_string("nonsense")


def test_sort_mode_from_string_round_trip():
    assert SortMode.from_string("edit_distance") == SortMode.EditDistance
    assert SortMode.from_string("coordinates") == SortMode.Coordinates


def test_sort_mode_invalid_token_raises():
    with pytest.raises(Exception):
        SortMode.from_string("nonsense")


def test_output_mode_from_string_round_trip():
    assert OutputMode.from_string("both") == OutputMode.Both
    with pytest.raises(Exception):
        OutputMode.from_string("nonsense")
