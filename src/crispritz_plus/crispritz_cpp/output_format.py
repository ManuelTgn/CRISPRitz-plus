"""Python wrapper over the C++ OutputFormat enum.

The members are the native C++ enum values, so they pass straight through
SearchConfiguration and the binding unchanged; the class adds friendly,
C++-validated string conversion on top.
"""

from crispritz_plus import _ternary_search_tree as tst  # type: ignore


class OutputFormat:
    """Serialization layout for the off-target table."""

    #: Canonical tab-separated values layout.
    Tsv = tst.OutputFormat.Tsv
    #: Legacy CRISPRitz "targets" column layout.
    Targets = tst.OutputFormat.Targets

    # Canonical lowercase tokens (mirror the C++ to_string()).
    _TOKENS = {Tsv: "tsv", Targets: "targets"}

    @classmethod
    def from_string(cls, name: str) -> "tst.OutputFormat":
        """Parse ``"tsv"`` / ``"targets"`` into the native enum value.

        Delegates to the C++ parser, so the accepted token set has a single
        source of truth.

        Raises
        ------
        TSTSearchError
            On an unrecognised token.
        """
        return tst.output_format_from_string(name)

    @classmethod
    def to_string(cls, value: "tst.OutputFormat") -> str:
        """Return the canonical lowercase token for a native enum value."""
        return cls._TOKENS[value]
