"""
Exception class for annotation-related errors in CRISPRitz-plus.

Defines :class:`CrispritzAnnotationError`, a thin subclass of
:class:`~crispritz_plus.crispritz_errors.CrispritzError` used to
distinguish failures that originate in the ``annotate-results`` pipeline
from other CRISPRitz-plus errors.

The separation lets callers catch annotation failures precisely without
accidentally masking unrelated :class:`CrispritzError` subtypes::

    try:
        annotate_results(...)
    except CrispritzAnnotationError as exc:
        # annotation-specific recovery
        ...
"""

from ..crispritz_errors import CrispritzError


class CrispritzAnnotationError(CrispritzError):
    """Exception raised for errors in the ``annotate-results`` pipeline.

    Inherits from :class:`~crispritz_plus.crispritz_errors.CrispritzError`
    and adds no new state or behaviour; its distinct type is the sole
    purpose, allowing callers to catch annotation failures independently of
    other CRISPRitz-plus error categories.

    Parameters
    ----------
    value : str
        Human-readable description of the error condition.

    Examples
    --------
    Raise and catch the exception::

        try:
            raise CrispritzAnnotationError("BED track is malformed")
        except CrispritzAnnotationError as exc:
            print(exc)   # BED track is malformed
    """

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
