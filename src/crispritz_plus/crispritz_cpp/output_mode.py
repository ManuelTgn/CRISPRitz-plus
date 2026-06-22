"""Python wrapper over the C++ OutputMode enum.

The members are the native C++ enum values, so they pass straight through
SearchConfiguration and the binding unchanged; the class adds friendly,
C++-validated string conversion on top.
"""

from crispritz_plus import _ternary_search_tree as tst  # type: ignore


class OutputMode:
    """Which result files a search run produces."""

    #: Write the off-target targets table only.
    TargetsOnly = tst.OutputMode.TargetsOnly
    #: Write the per-guide profile files only.
    ProfileOnly = tst.OutputMode.ProfileOnly
    #: Write both the targets table and the profile files.
    Both = tst.OutputMode.Both

    # Canonical lowercase tokens (mirror the C++ to_string()).
    _TOKENS = {TargetsOnly: "targets", ProfileOnly: "profile", Both: "both"}

    @classmethod
    def from_string(cls, name: str) -> "tst.OutputMode":
        """Parse ``"targets"`` / ``"profile"`` / ``"both"`` into the native enum.

        Delegates to the C++ parser, so the accepted token set has a single
        source of truth.

        Raises
        ------
        TSTSearchError
            On an unrecognised token.
        """
        return tst.output_mode_from_string(name)

    @classmethod
    def to_string(cls, value: "tst.OutputMode") -> str:
        """Return the canonical lowercase token for a native enum value."""
        return cls._TOKENS[value]
