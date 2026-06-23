""" """

from crispritz_plus import _ternary_search_tree as tst  # type: ignore


class BulgeMode:
    """Whether one off-target alignment may combine DNA and RNA bulges."""

    #: A hit may contain both DNA and RNA bulges at once, each up to its
    #: per-kind budget. The default; the legacy token ``"both"`` also parses to
    #: this value.
    MixedBulges = tst.BulgeMode.MixedBulges
    #: A hit may use only one bulge kind — DNA or RNA, never both in the same
    #: target. Mismatches combine freely with either.
    SingleBulgeType = tst.BulgeMode.SingleBulgeType

    # Canonical lowercase tokens (mirror the C++ parser's primary spellings).
    _TOKENS = {MixedBulges: "mixed", SingleBulgeType: "single"}

    @classmethod
    def from_string(cls, name: str) -> "tst.BulgeMode":
        """Parse ``"mixed"`` / ``"single"`` (or legacy ``"both"``) into the native enum.

        Delegates to the C++ parser, so the accepted token set has a single
        source of truth.

        Raises
        ------
        TSTSearchError
            On an unrecognised token.
        """
        return tst.bulge_mode_from_string(name)

    @classmethod
    def to_string(cls, value: "tst.BulgeMode") -> str:
        """Return the canonical lowercase token for a native enum value."""
        return cls._TOKENS[value]
