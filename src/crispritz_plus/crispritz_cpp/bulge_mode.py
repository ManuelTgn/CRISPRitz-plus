"""Python wrapper over the C++ ``BulgeMode`` enum.

Controls whether a single off-target alignment may mix DNA and RNA bulges.
The members are the native C++ enum values, so they pass straight through the
search entry point unchanged; the facade adds friendly, C++-validated string
conversion on top.
"""

from crispritz_plus import _ternary_search_tree as tst  # type: ignore


class BulgeMode:
    """Whether one off-target alignment may combine DNA and RNA bulges.

    The two members mirror the native ``tst.BulgeMode`` values.  Mismatches
    combine freely with either bulge kind regardless of the mode.
    """

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
        """Parse ``"mixed"`` / ``"single"`` (or the legacy ``"both"``) into the native enum value.

        Delegates to the C++ parser, so the accepted token set has a single
        source of truth.

        Parameters
        ----------
        name : str
            The token to parse.

        Returns
        -------
        "tst.BulgeMode"
            The corresponding native C++ enum value.

        Raises
        ------
        TSTSearchError
            On an unrecognised token.
        """
        return tst.bulge_mode_from_string(name)

    @classmethod
    def to_string(cls, value: "tst.BulgeMode") -> str:
        """Return the canonical lowercase token for a native enum value.

        Parameters
        ----------
        value : "tst.BulgeMode"
            A native C++ enum value.

        Returns
        -------
        str
            The canonical lowercase token for *value*.
        """
        return cls._TOKENS[value]
