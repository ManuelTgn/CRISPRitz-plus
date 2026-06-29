"""
FASTA / VCF input pairing for the CRISPRitz-plus enrichment pipeline.

Defines :class:`EnrichPair`, a small validated container that holds the two
inputs required to enrich a reference genome with variants: a FASTA file and
a VCF file.  Each path is validated on assignment so that downstream
enrichment code can assume well-formed inputs.
"""

import os


from ..exception_handlers import exception_handler

from .crispritz_enrichment_error import EnrichmentPairError


class EnrichPair:
    """Validated FASTA / VCF input pair for one enrichment unit.

    Holds the two file paths consumed by the variant-enrichment step.  Both
    fields are exposed as properties whose setters validate the supplied
    value, raising :class:`~.crispritz_enrichment_error.EnrichmentPairError`
    (via :func:`~crispritz_plus.exception_handlers.exception_handler`) on a
    type or emptiness violation.

    Parameters
    ----------
    debug : bool
        When *True*, validation errors propagate with a full traceback
        instead of a formatted user-facing message.

    Attributes
    ----------
    fasta : str
        Path to the reference FASTA file (validated non-empty ``str``).
    vcf : str
        Path to the VCF file (validated ``str``).
    """

    def __init__(self, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._fasta: str = ""
        self._vcf: str = ""

    def __repr__(self) -> str:
        """Return a detailed string representation for debugging."""
        return f"EnrichPair(fasta={self._fasta!r}, vcf={self._vcf!r})"

    def __str__(self) -> str:
        """Return a human-readable string representation."""
        fasta_str = self._fasta if self._fasta is not None else "not set"
        vcf_str = self._vcf if self._vcf is not None else "not set"
        return f"EnrichPair: FASTA={fasta_str}, VCF={vcf_str}"

    @property
    def fasta(self) -> str:
        """str: Path to the reference FASTA file."""
        return self._fasta

    @fasta.setter
    def fasta(self, value: str) -> None:
        """Set the FASTA path, validating that it is a non-empty string.

        Parameters
        ----------
        value : str
            Path to the reference FASTA file.

        Raises
        ------
        EnrichmentPairError
            If *value* is not a ``str`` or is an empty string.
        """
        if not isinstance(value, str) or not value:
            exception_handler(
                EnrichmentPairError,
                "FASTA file must be a non-empty str, got "
                f"{type(value).__name__} instead",
                os.EX_DATAERR,
                self._debug,
            )
        self._fasta = value

    @property
    def vcf(self) -> str:
        """str: Path to the VCF file."""
        return self._vcf

    @vcf.setter
    def vcf(self, value: str) -> None:
        """Set the VCF path, validating that it is a string.

        Unlike :attr:`fasta`, an empty string is accepted here; only the
        type is enforced.

        Parameters
        ----------
        value : str
            Path to the VCF file.

        Raises
        ------
        EnrichmentPairError
            If *value* is not a ``str``.
        """
        if not isinstance(value, str):
            exception_handler(
                EnrichmentPairError,
                "VCF file must be a non-empty str, got "
                f"{type(value).__name__} instead",
                os.EX_DATAERR,
                self._debug,
            )
        self._vcf = value
