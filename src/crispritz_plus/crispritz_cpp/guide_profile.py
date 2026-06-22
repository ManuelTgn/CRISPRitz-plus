"""Python wrapper over the C++ GuideProfile.

A GuideProfile holds a guide's accumulated profiling statistics; it is produced
C++-side and consumed C++-side (it is passed back to the profile-merge entry
point as an opaque value). This wrapper takes the native object, delegates the
read-only fields exposed across the binding, and exposes the underlying object
via :pyattr:`native` so it can be handed back across the boundary.

Only the fields the binding surfaces (guide identity and on-target counts) are
delegated; the full counter matrices stay C++-internal and are written by the
ProfileWriter, never inspected from Python.
"""

from crispritz_plus import _ternary_search_tree as tst  # type: ignore


class GuideProfile:
    """A single guide's profiling statistics (post-accumulation)."""

    __slots__ = ("_profile",)

    def __init__(self, profile: object) -> None:
        # Accept a native C++ GuideProfile or another GuideProfile wrapper.
        self._profile = getattr(profile, "native", profile)

    # ------------------------------------------------------------------
    # Access to the wrapped C++ object
    # ------------------------------------------------------------------

    @property
    def native(self) -> "tst.GuideProfile":
        """The underlying C++ ``GuideProfile`` (for the binding boundary)."""
        return self._profile

    # ------------------------------------------------------------------
    # Delegating read-only accessors
    # ------------------------------------------------------------------

    @property
    def guide(self) -> str:
        """Guide body sequence (PAM bases excluded)."""
        return self._profile.guide  # type: ignore

    @property
    def guide_len(self) -> int:
        """Length of the guide body."""
        return self._profile.guide_len  # type: ignore

    @property
    def ont_count(self) -> int:
        """On-target hit count (mismatch-only channel: 0 mismatches, 0 bulge)."""
        return self._profile.ont_count  # type: ignore

    @property
    def ont_count_complete(self) -> int:
        """On-target hit count across all channels (any 0-mismatch hit)."""
        return self._profile.ont_count_complete  # type: ignore

    def __repr__(self) -> str:
        return (
            "GuideProfile("
            f"guide={self._profile.guide!r}, "  # type: ignore
            f"guide_len={self._profile.guide_len}, "  # type: ignore
            f"ont_count={self._profile.ont_count}, "  # type: ignore
            f"ont_count_complete={self._profile.ont_count_complete})"  # type: ignore
        )
