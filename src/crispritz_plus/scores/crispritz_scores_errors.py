"""Exception types for the CRISPRitz-plus scoring stage.

Defines :class:`CfdScoreError`, a subclass of
:class:`~crispritz_plus.crispritz_errors.CrispritzError` raised when the CFD
model files cannot be loaded or a scoring operation fails, so scoring failures
can be caught independently of other CRISPRitz-plus errors.
"""

from ..crispritz_errors import CrispritzError


class CfdScoreError(CrispritzError):
    """Exception raised for CFD scoring failures.

    Inherits from :class:`~crispritz_plus.crispritz_errors.CrispritzError` and
    adds no new behaviour; its distinct type lets callers catch CFD scoring
    failures (e.g. missing model pickles) independently.

    Parameters
    ----------
    value : str
        Human-readable description of the error condition.
    """

    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
