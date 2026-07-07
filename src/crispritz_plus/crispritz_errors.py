"""Exception hierarchy for CRISPRitz-plus.

Defines :class:`CrispritzError`, the root of every custom exception in the
package, together with the stage- and component-specific subclasses raised by
:func:`~crispritz_plus.exception_handlers.exception_handler`.  Grouping errors
under one root lets callers catch the whole family with a single
``except CrispritzError`` while still allowing fine-grained handling of
individual subtypes.

Hierarchy
---------
CrispritzError
    Root of the hierarchy.
    ├── GenomeIoError
    │   ├── GenomeReaderError
    │   └── GenomeWriterError
    ├── CrispritzPamError
    ├── CrispritzTstError
    ├── CrispritzSearchError
    └── CrispritzGuideError
"""


class CrispritzError(Exception):
    """Root exception type for all CRISPRitz-plus errors.

    Parameters
    ----------
    value : str
        Human-readable description of the error condition.
    """

    def __init__(self, value: str):
        # initialize exception object when raised
        self._value = value  # error message or error related info

    def __str__(self):
        return repr(self._value)  # string representation for the exception


class GenomeIoError(CrispritzError):
    """Base error for genome FASTA input/output failures.

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


class GenomeReaderError(GenomeIoError):
    """Error raised while reading a genome FASTA file.

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


class GenomeWriterError(GenomeIoError):
    """Error raised while writing a genome FASTA file.

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


class CrispritzPamError(CrispritzError):
    """Error raised while parsing or validating a PAM.

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


class CrispritzTstError(CrispritzError):
    """Error raised during Ternary Search Tree operations.

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


class CrispritzSearchError(CrispritzError):
    """Error raised during the off-target search stage.

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


class CrispritzGuideError(CrispritzError):
    """Error raised while parsing or validating guides.

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
