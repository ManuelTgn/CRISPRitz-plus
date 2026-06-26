"""
Exception class for report-generation errors in CRISPRitz-plus.

Defines :class:`CrispritzReportError`, a thin subclass of
:class:`~crispritz_plus.crispritz_errors.CrispritzError` used to
distinguish failures that originate in the ``generate-report`` pipeline
from other CRISPRitz-plus errors.

The separation lets callers catch report failures precisely without
accidentally masking unrelated :class:`CrispritzError` subtypes::

    try:
        generate_report_cli(args)
    except CrispritzReportError as exc:
        # report-specific recovery
        ...
"""

from ..crispritz_errors import CrispritzError


class CrispritzReportError(CrispritzError):
    """Exception raised for errors in the ``generate-report`` pipeline.

    Inherits from :class:`~crispritz_plus.crispritz_errors.CrispritzError`
    and adds no new state or behaviour; its distinct type is the sole
    purpose, allowing callers to catch report-generation failures
    independently of other CRISPRitz-plus error categories.

    Parameters
    ----------
    value : str
        Human-readable description of the error condition.

    Examples
    --------
    Raise and catch the exception::

        try:
            raise CrispritzReportError("annotated TSV is empty")
        except CrispritzReportError as exc:
            print(exc)   # annotated TSV is empty
    """

    def __init__(self, value: str):
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
