"""
tests/core/test_variants.py
============================
Unit tests for ``crispritz_plus.enrichment.variants``.

Covers the value objects (``IndelInfo``, ``IndelPair``) and the carrier
types (``Snp``/``Snps``, ``Indel``/``Indels``, ``IndelsSet``).
"""

from __future__ import annotations

import pytest

from crispritz_plus.variants import (
    IndelInfo,
    IndelPair,
    Snp,
    Snps,
    Indel,
    Indels,
    IndelsSet,
)


# ===========================================================================
# IndelInfo / IndelPair (dataclasses)
# ===========================================================================


class TestIndelInfo:
    def test_field_access(self):
        info = IndelInfo(idx=1, start=10, stop=20)
        assert info.idx == 1
        assert info.start == 10
        assert info.stop == 20

    def test_equality(self):
        assert IndelInfo(1, 10, 20) == IndelInfo(1, 10, 20)

    def test_inequality(self):
        assert IndelInfo(1, 10, 20) != IndelInfo(2, 10, 20)


class TestIndelPair:
    def test_field_access(self):
        pair = IndelPair(refseq=list("ACGT"), indelseq=list("ACXGT"))
        assert pair.refseq == list("ACGT")
        assert pair.indelseq == list("ACXGT")

    def test_equality(self):
        a = IndelPair(list("AC"), list("AC"))
        b = IndelPair(list("AC"), list("AC"))
        assert a == b


# ===========================================================================
# Snp
# ===========================================================================


class TestSnp:
    def test_properties(self):
        s = Snp(pos=100, ref="A", alt="G", gtidx=1)
        assert s.pos == 100
        assert s.ref == "A"
        assert s.alt == "G"
        assert s.gtidx == 1

    def test_repr_contains_fields(self):
        s = Snp(pos=5, ref="C", alt="T", gtidx=2)
        r = repr(s)
        assert "pos=5" in r
        assert "'C'" in r
        assert "'T'" in r


# ===========================================================================
# Snps (collection)
# ===========================================================================


class TestSnps:
    def test_empty_collection_is_falsy(self):
        assert bool(Snps()) is False

    def test_nonempty_collection_is_truthy(self):
        snps = Snps([Snp(0, "A", "G", 1)])
        assert bool(snps) is True

    def test_len(self):
        snps = Snps([Snp(0, "A", "G", 1), Snp(0, "A", "T", 2)])
        assert len(snps) == 2

    def test_iteration_yields_snp_objects(self):
        items = [Snp(0, "A", "G", 1), Snp(0, "A", "T", 2)]
        snps = Snps(items)
        assert list(snps) == items

    def test_default_constructor_is_empty(self):
        assert len(Snps()) == 0

    def test_add_appends_snp(self):
        snps = Snps()
        snps.add(Snp(10, "C", "G", 1))
        assert len(snps) == 1

    def test_add_rejects_non_snp(self):
        snps = Snps()
        with pytest.raises(TypeError):
            snps.add("not a snp")

    def test_items_returns_copy_not_reference(self):
        snps = Snps([Snp(0, "A", "G", 1)])
        copy = snps.items
        copy.append(Snp(99, "T", "A", 2))
        assert len(snps) == 1  # original collection unaffected

    def test_alts_in_order(self):
        snps = Snps([Snp(0, "A", "G", 1), Snp(0, "A", "T", 2)])
        assert snps.alts == ["G", "T"]

    def test_gtidxs_in_order(self):
        snps = Snps([Snp(0, "A", "G", 1), Snp(0, "A", "T", 2)])
        assert snps.gtidxs == [1, 2]

    def test_pos_returns_shared_position(self):
        snps = Snps([Snp(50, "A", "G", 1), Snp(50, "A", "T", 2)])
        assert snps.pos == 50

    def test_pos_raises_on_empty(self):
        with pytest.raises(AssertionError):
            _ = Snps().pos

    def test_ref_returns_shared_reference(self):
        snps = Snps([Snp(50, "A", "G", 1), Snp(50, "A", "T", 2)])
        assert snps.ref == "A"

    def test_ref_raises_on_empty(self):
        with pytest.raises(AssertionError):
            _ = Snps().ref

    def test_repr_contains_count(self):
        snps = Snps([Snp(0, "A", "G", 1)])
        assert "n=1" in repr(snps)


# ===========================================================================
# Indel
# ===========================================================================


class TestIndel:
    def test_properties(self):
        i = Indel(pos=200, ref="AT", alt="A", gtidx=1)
        assert i.pos == 200
        assert i.ref == "AT"
        assert i.alt == "A"
        assert i.gtidx == 1

    def test_repr_contains_fields(self):
        i = Indel(pos=7, ref="A", alt="ATG", gtidx=3)
        r = repr(i)
        assert "pos=7" in r
        assert "'ATG'" in r

    def test_is_symbolic_false_for_literal_alleles(self):
        i = Indel(pos=0, ref="A", alt="ATG", gtidx=1)
        assert i.is_symbolic() is False

    def test_is_symbolic_true_when_alt_is_symbolic(self):
        i = Indel(pos=0, ref="A", alt="<DEL>", gtidx=1)
        assert i.is_symbolic() is True

    def test_is_symbolic_true_when_ref_is_symbolic(self):
        i = Indel(pos=0, ref="<INV>", alt="A", gtidx=1)
        assert i.is_symbolic() is True

    def test_is_symbolic_true_when_both_symbolic(self):
        i = Indel(pos=0, ref="<INV>", alt="<DEL>", gtidx=1)
        assert i.is_symbolic() is True


# ===========================================================================
# Indels (collection)
# ===========================================================================


class TestIndels:
    def test_empty_is_falsy(self):
        assert bool(Indels()) is False

    def test_nonempty_is_truthy(self):
        assert bool(Indels([Indel(0, "A", "AT", 1)])) is True

    def test_len(self):
        indels = Indels([Indel(0, "A", "AT", 1), Indel(0, "A", "AG", 2)])
        assert len(indels) == 2

    def test_iteration_yields_indel_objects(self):
        items = [Indel(0, "A", "AT", 1)]
        indels = Indels(items)
        assert list(indels) == items

    def test_add_appends_indel(self):
        indels = Indels()
        indels.add(Indel(0, "A", "AT", 1))
        assert len(indels) == 1

    def test_add_rejects_non_indel(self):
        indels = Indels()
        with pytest.raises(TypeError):
            indels.add(("A", "AT"))

    def test_items_returns_copy(self):
        indels = Indels([Indel(0, "A", "AT", 1)])
        copy = indels.items
        copy.append(Indel(99, "G", "GA", 2))
        assert len(indels) == 1

    def test_alts_in_order(self):
        indels = Indels([Indel(0, "A", "AT", 1), Indel(0, "A", "AGG", 2)])
        assert indels.alts() == ["AT", "AGG"]

    def test_gtidxs_in_order(self):
        indels = Indels([Indel(0, "A", "AT", 1), Indel(0, "A", "AGG", 2)])
        assert indels.gtidxs() == [1, 2]

    def test_pos_returns_shared_position(self):
        indels = Indels([Indel(75, "A", "AT", 1), Indel(75, "A", "AGG", 2)])
        assert indels.pos() == 75

    def test_pos_raises_on_empty(self):
        with pytest.raises(AssertionError):
            Indels().pos()

    def test_ref_returns_shared_reference(self):
        indels = Indels([Indel(75, "GC", "G", 1)])
        assert indels.ref() == "GC"

    def test_ref_raises_on_empty(self):
        with pytest.raises(AssertionError):
            Indels().ref()

    def test_non_symbolic_filters_symbolic_alleles(self):
        indels = Indels(
            [
                Indel(0, "A", "AT", 1),  # literal
                Indel(0, "A", "<DEL>", 2),  # symbolic
            ]
        )
        filtered = indels.non_symbolic()
        assert len(filtered) == 1
        assert filtered.items[0].alt == "AT"

    def test_non_symbolic_returns_new_indels_instance(self):
        indels = Indels([Indel(0, "A", "AT", 1)])
        filtered = indels.non_symbolic()
        assert isinstance(filtered, Indels)
        assert filtered is not indels

    def test_non_symbolic_all_literal_keeps_all(self):
        indels = Indels([Indel(0, "A", "AT", 1), Indel(0, "A", "AGG", 2)])
        assert len(indels.non_symbolic()) == 2

    def test_non_symbolic_all_symbolic_returns_empty(self):
        indels = Indels([Indel(0, "A", "<DEL>", 1), Indel(0, "<INV>", "A", 2)])
        assert len(indels.non_symbolic()) == 0

    def test_repr_contains_count(self):
        indels = Indels([Indel(0, "A", "AT", 1)])
        assert "n=1" in repr(indels)


# ===========================================================================
# IndelsSet
# ===========================================================================


class TestIndelsSet:
    def test_first_push_starts_at_zero(self):
        s = IndelsSet(debug=True)
        info = s.push(list("ACGT"))
        assert info.start == 0
        assert info.stop == 4

    def test_second_push_leaves_one_base_gap(self):
        s = IndelsSet(debug=True)
        s.push(list("ACGT"))  # [0, 4)
        info2 = s.push(list("TT"))  # should start at 5 (gap)
        assert info2.start == 5
        assert info2.stop == 7

    def test_idx_starts_at_two_after_first_push(self):
        # self._i starts at 1, incremented to 2 before returning
        s = IndelsSet(debug=True)
        info = s.push(list("AC"))
        assert info.idx == 2

    def test_idx_increments_each_push(self):
        s = IndelsSet(debug=True)
        info1 = s.push(list("A"))
        info2 = s.push(list("C"))
        info3 = s.push(list("G"))
        assert [info1.idx, info2.idx, info3.idx] == [2, 3, 4]

    def test_sequences_accumulate_in_order(self):
        s = IndelsSet(debug=True)
        s.push(list("AC"))
        s.push(list("GT"))
        assert s.sequences == [list("AC"), list("GT")]

    def test_start_i_reflects_cursor_after_pushes(self):
        s = IndelsSet(debug=True)
        s.push(list("AC"))  # cursor -> 2+1=3
        s.push(list("G"))  # cursor -> 3+1+1=5
        assert s.start_i == 5

    def test_empty_sequence_push_still_advances_by_gap(self):
        s = IndelsSet(debug=True)
        info = s.push([])
        assert info.start == 0
        assert info.stop == 0
        assert s.start_i == 1  # one-base gap even for empty sequence

    def test_returns_indelinfo_instance(self):
        s = IndelsSet(debug=True)
        assert isinstance(s.push(list("A")), IndelInfo)

    def test_three_consecutive_pushes_no_overlap(self):
        s = IndelsSet(debug=True)
        spans = [s.push(list("A" * n)) for n in (3, 5, 2)]
        # No span should overlap with the next: stop + 1 == next start
        for cur, nxt in zip(spans, spans[1:]):
            assert nxt.start == cur.stop + 1
