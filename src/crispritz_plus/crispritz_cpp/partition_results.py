"""Python wrapper over the C++ ``PartitionResult``.

A ``PartitionResult`` is produced by the C++ search executor (one per
partition); Python never builds one from scratch.  This wrapper takes the
native object and delegates its read-only fields, exposing the underlying
object via :attr:`native` for the rare places that must hand the real type
back across the pybind11 boundary (e.g. the profile-merge entry point, which
consumes the native ``profiles``).
"""

from typing import List


from crispritz_plus import _ternary_search_tree as tst  # type: ignore

from .guide_profile import GuideProfile


class PartitionResult:
    """Outcome of searching a single partition: counts plus per-guide profiles.

    Wraps a native ``tst.PartitionResult`` and delegates its read-only fields.

    Parameters
    ----------
    result : object
        A native C++ ``PartitionResult`` or another ``PartitionResult``
        wrapper; construction is idempotent (see :meth:`__init__`).
    """

    __slots__ = ("_result",)

    def __init__(self, result: object) -> None:
        # Accept a native C++ PartitionResult or another PartitionResult wrapper.
        self._result = getattr(result, "native", result)

    # ==========================================================================
    # Access to the wrapped C++ object
    # ==========================================================================

    @property
    def native(self) -> "tst.PartitionResult":
        """The underlying C++ ``PartitionResult`` (for the binding boundary)."""
        return self._result

    # ==========================================================================
    # Delegating read-only accessors
    # ==========================================================================

    @property
    def source_path(self) -> str:
        """Path of the ``.bin`` partition these results came from."""
        return self._result.source_path  # type: ignore

    @property
    def shard_path(self) -> str:
        """Shard file the targets were written to ("" when targets disabled)."""
        return self._result.shard_path  # type: ignore

    @property
    def total_hits(self) -> int:
        """Total hits found across all guides in this partition."""
        return self._result.total_hits  # type: ignore

    @property
    def rows_written(self) -> int:
        """Rows written to the shard (== total_hits when targets enabled)."""
        return self._result.rows_written  # type: ignore

    @property
    def profiles(self) -> List[GuideProfile]:
        """One native ``GuideProfile`` per guide (empty when profiling disabled)."""
        return self._result.profiles  # type: ignore

    def __repr__(self) -> str:
        return (
            "PartitionResult("
            f"source_path={self._result.source_path!r}, "  # type: ignore
            f"shard_path={self._result.shard_path!r}, "  # type: ignore
            f"total_hits={self._result.total_hits}, "  # type: ignore
            f"rows_written={self._result.rows_written}, "  # type: ignore
            f"profiles={len(self._result.profiles)})"  # type: ignore
        )
