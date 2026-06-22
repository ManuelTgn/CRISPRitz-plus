"""Python wrapper over the C++ SortMode enum.

The members are the native C++ enum values, so they pass straight through the
merge entry point unchanged; the class adds friendly, C++-validated string
conversion on top. SortMode is API-only (it is not a CLI flag).
"""

from crispritz_plus import _ternary_search_tree as tst  # type: ignore


class SortMode:
    """Ordering applied to the final off-target table."""

    #: edits ascending, then mismatches, then bulges, then CFD descending
    #: (NA last). The default.
    EditDistance = tst.SortMode.EditDistance
    #: contig (lexicographic) ascending, then position ascending.
    Coordinates = tst.SortMode.Coordinates

    # Canonical lowercase tokens (mirror the C++ to_string()).
    _TOKENS = {EditDistance: "edit_distance", Coordinates: "coordinates"}

    @classmethod
    def from_string(cls, name: str) -> "tst.SortMode":
        """Parse ``"edit_distance"`` / ``"coordinates"`` into the native enum.

        Delegates to the C++ parser, so the accepted token set has a single
        source of truth.

        Raises
        ------
        TSTSearchError
            On an unrecognised token.
        """
        return tst.sort_mode_from_string(name)

    @classmethod
    def to_string(cls, value: "tst.SortMode") -> str:
        """Return the canonical lowercase token for a native enum value."""
        return cls._TOKENS[value]
