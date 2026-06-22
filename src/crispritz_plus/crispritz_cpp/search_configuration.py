"""High-level Python wrapper around the C++ SearchConfiguration.

The low-level, validated configuration object lives in the C++ extension. This
wrapper provides a friendlier Python surface over it:

  * its constructor accepts the lowercase CLI tokens for ``output_format`` and
    ``output_mode`` (e.g. ``"tsv"``, ``"both"``) as well as the bound enum
    values, parsing tokens through the C++ ``*_from_string`` helpers so the
    canonical token set stays defined in one place (C++);
  * it exposes the underlying C++ object via :pyattr:`native` for the rare
    places that must hand the real type across the pybind11 boundary
    (the search-executor entry point).

All validation still happens in C++ at construction time; the wrapper holds no
state of its own beyond the wrapped object and delegates every accessor to it.
"""

from typing import Union

from crispritz_plus import _ternary_search_tree as tst  # type: ignore

from .output_format import OutputFormat
from .output_mode import OutputMode


class SearchConfiguration:

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

    # ------------------------------------------------------------------
    # Access to the wrapped C++ object
    # ------------------------------------------------------------------

    @property
    def native(self) -> "tst.SearchConfiguration":
        """The underlying C++ ``SearchConfiguration`` (for the binding boundary)."""
        return self._config

    # ------------------------------------------------------------------
    # Delegating read-only accessors
    # ------------------------------------------------------------------

    @property
    def max_mismatches(self) -> int:
        return self._config.max_mismatches

    @property
    def max_bulges_dna(self) -> int:
        return self._config.max_bulges_dna

    @property
    def max_bulges_rna(self) -> int:
        return self._config.max_bulges_rna

    @property
    def threads(self) -> int:
        return self._config.threads

    @property
    def output_format(self) -> str:
        return OutputFormat.to_string(self._config.output_format)

    @property
    def output_mode(self) -> str:
        return OutputMode.to_string(self._config.output_format)

    @property
    def max_total_edits(self) -> int:
        return self._config.max_total_edits

    @property
    def write_targets(self) -> bool:
        return self._config.write_targets

    @property
    def write_profile(self) -> bool:
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
