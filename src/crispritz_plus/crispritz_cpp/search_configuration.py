"""High-level Python wrapper around the C++ ``SearchConfiguration``.

The low-level, validated configuration object lives in the C++ extension.
This wrapper provides a friendlier Python surface over it:

* its constructor accepts the lowercase CLI tokens for ``output_format`` and
  ``output_mode`` (e.g. ``"tsv"``, ``"both"``) as well as the bound enum
  values, parsing tokens through the C++ ``*_from_string`` helpers so the
  canonical token set stays defined in one place (C++);
* it exposes the underlying C++ object via :attr:`native` for the rare places
  that must hand the real type across the pybind11 boundary (the
  search-executor entry point).

All validation still happens in C++ at construction time; the wrapper holds no
state of its own beyond the wrapped object and delegates every accessor to it.
"""

from typing import Union


from crispritz_plus import _ternary_search_tree as tst  # type: ignore

from .output_format import OutputFormat
from .output_mode import OutputMode


class SearchConfiguration:
    """Validated parameters for one off-target search run.

    Thin facade over the native ``tst.SearchConfiguration``.  String tokens
    for *output_format* and *output_mode* are resolved to native enum values
    at construction; every accessor delegates to the wrapped C++ object.

    Parameters
    ----------
    max_mismatches : int
        Maximum number of mismatches permitted in a hit.
    max_bulges_dna : int
        Maximum number of DNA bulges permitted in a hit.
    max_bulges_rna : int
        Maximum number of RNA bulges permitted in a hit.
    threads : int
        Number of threads the C++ search may use.
    output_format : Union[str, OutputFormat], optional
        Off-target table layout, as a token (e.g. ``"tsv"``) or an
        :class:`OutputFormat` value. Defaults to ``"tsv"``.
    output_mode : Union[str, OutputMode], optional
        Which result files to produce, as a token (e.g. ``"both"``) or an
        :class:`OutputMode` value. Defaults to ``"both"``.
    """

    __slots__ = ("_config",)

    def __init__(
        self,
        max_mismatches: int,
        max_bulges_dna: int,
        max_bulges_rna: int,
        threads: int,
        output_format: Union[str, OutputFormat] = "tsv",
        output_mode: Union[str, OutputMode] = "both",
    ) -> None:
        fmt = (
            OutputFormat.from_string(output_format)
            if isinstance(output_format, str)
            else output_format
        )
        mode = (
            OutputMode.from_string(output_mode)
            if isinstance(output_mode, str)
            else output_mode
        )
        self._config = tst.SearchConfiguration(
            max_mismatches, max_bulges_dna, max_bulges_rna, threads, fmt, mode
        )

    # ==========================================================================
    # Access to the wrapped C++ object
    # ==========================================================================

    @property
    def native(self) -> "tst.SearchConfiguration":
        """The underlying C++ ``SearchConfiguration`` (for the binding boundary)."""
        return self._config

    # ==========================================================================
    # Delegating read-only accessors
    # ==========================================================================

    @property
    def max_mismatches(self) -> int:
        """int: Maximum number of mismatches permitted in a hit."""
        return self._config.max_mismatches

    @property
    def max_bulges_dna(self) -> int:
        """int: Maximum number of DNA bulges permitted in a hit."""
        return self._config.max_bulges_dna

    @property
    def max_bulges_rna(self) -> int:
        """int: Maximum number of RNA bulges permitted in a hit."""
        return self._config.max_bulges_rna

    @property
    def threads(self) -> int:
        """int: Number of threads the C++ search may use."""
        return self._config.threads

    @property
    def output_format(self) -> str:
        """str: Off-target table layout as its canonical lowercase token."""
        return OutputFormat.to_string(self._config.output_format)

    @property
    def output_mode(self) -> str:
        """str: Which result files are produced, as a canonical lowercase token."""
        return OutputMode.to_string(self._config.output_format)

    @property
    def max_total_edits(self) -> int:
        """int: Maximum combined edits (mismatches + bulges) permitted in a hit."""
        return self._config.max_total_edits

    @property
    def write_targets(self) -> bool:
        """bool: Whether the targets table is written under the selected output mode."""
        return self._config.write_targets

    @property
    def write_profile(self) -> bool:
        """bool: Whether per-guide profiles are written under the selected output mode."""
        return self._config.write_profile

    def __repr__(self) -> str:
        return (
            "SearchConfiguration("
            f"max_mismatches={self.max_mismatches}, "
            f"max_bulges_dna={self.max_bulges_dna}, "
            f"max_bulges_rna={self.max_bulges_rna}, "
            f"threads={self.threads}, "
            f"write_targets={self.write_targets}, "
            f"write_profile={self.write_profile})"
        )
