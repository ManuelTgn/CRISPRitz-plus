"""Python wrapper over the C++ ``OutputMode`` enum.

The members are the native C++ enum values, so they pass straight through
:class:`~crispritz_plus.crispritz_cpp.search_configuration.SearchConfiguration`
and the binding unchanged; the facade adds friendly, C++-validated string
conversion on top.
"""

from crispritz_plus import _ternary_search_tree as tst  # type: ignore


class OutputMode:
    """Which result files a search run produces.

    The three members mirror the native ``tst.OutputMode`` values: the targets
    table only, the per-guide profiles only, or both.
    """

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
        """Parse ``"targets"`` / ``"profile"`` / ``"both"`` into the native enum value.

        Delegates to the C++ parser, so the accepted token set has a single
        source of truth.

        Parameters
        ----------
        name : str
            The token to parse.

        Returns
        -------
        "tst.OutputMode"
            The corresponding native C++ enum value.

        Raises
        ------
        TSTSearchError
            On an unrecognised token.
        """
        return tst.output_mode_from_string(name)

    @classmethod
    def to_string(cls, value: "tst.OutputMode") -> str:
        """Return the canonical lowercase token for a native enum value.

        Parameters
        ----------
        value : "tst.OutputMode"
            A native C++ enum value.

        Returns
        -------
        str
            The canonical lowercase token for *value*.
        """
        return cls._TOKENS[value]
