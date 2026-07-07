"""PAM parsing, geometry, and Cas-system recognition for CRISPRitz-plus.

Defines :class:`PAM`, which reads a PAM specification from file, derives its
length, orientation (upstream vs. downstream of the spacer), and implied guide
length, and matches it against the known PAM catalogues to identify the Cas
system.

Module-level constants
----------------------
CASXPAM, CPF1PAM, SACAS9PAM, SPCAS9PAM, XCAS9PAM : List[str]
    Known PAM motifs for each supported Cas system.
CASX, CPF1, SACAS9, SPCAS9, XCAS9 : int
    Integer identifiers for the corresponding Cas systems; ``-1`` denotes an
    unrecognised system.
"""

from typing import Tuple

import os


from .crispritz_errors import CrispritzPamError
from .exception_handlers import exception_handler


# ==============================================================================
# PAM domain constants
# ==============================================================================

#: Known PAM motif(s) for the CasX (Cas12e) system.
CASXPAM = ["TTCN"]

#: Known PAM motifs for the Cpf1 (Cas12a) system (upstream / 5' PAMs).
CPF1PAM = [
    "TTN",
    "TTTN",
    "TYCV",
    "TATV",
    "TTTV",
    "TTTR",
    "ATTN",
    "TTTA",
    "TCTA",
    "TCCA",
    "CCCA",
    "YTTV",
    "TTYN",
]

#: Known PAM motifs for the SaCas9 system.
SACAS9PAM = ["NNGRRT", "NNNRRT"]

#: Known PAM motifs for the SpCas9 system.
SPCAS9PAM = ["NGG", "NGA", "NRG", "NGC"]

#: Known PAM motifs for the xCas9 system.
XCAS9PAM = ["NGK", "NGN", "NNG"]

#: Integer identifier for the CasX (Cas12e) system.
CASX = 0

#: Integer identifier for the Cpf1 (Cas12a) system.
CPF1 = 1

#: Integer identifier for the SaCas9 system.
SACAS9 = 2

#: Integer identifier for the SpCas9 system.
SPCAS9 = 3

#: Integer identifier for the xCas9 system.
XCAS9 = 4


# ==============================================================================
# PAM class definition
# ==============================================================================


class PAM:
    """A parsed PAM with orientation, size, and Cas-system classification.

    Reads the PAM sequence and a signed size from file (a negative size marks
    an upstream PAM), derives the spacer length, isolates the PAM motif, and
    classifies it against the known Cas-system catalogues.

    Parameters
    ----------
    fname : str
        Path to the PAM specification file.
    debug : bool
        When *True*, parsing errors propagate with a full traceback.
    """

    def __init__(self, fname: str, debug: bool) -> None:
        self._debug = debug  # store debug flag
        # set pam sequence and size
        pam_seq, size = self._read_pam_file(fname)
        self._upstream = size < 0  # pam is upstream
        self._size = abs(size)  # adjust size
        self._guide_size = len(pam_seq) - self._size
        self._pam_seq = self._refine_pam(pam_seq)
        self._assess_cas_system()

    def __len__(self) -> int:
        return self._size

    def _read_pam_file(self, fname: str) -> Tuple[str, int]:
        """Read the PAM sequence and signed size from file.

        The file's first line holds the PAM sequence and a size; a negative
        size indicates the PAM lies upstream of the guide.

        Parameters
        ----------
        fname : str
            Path to the PAM specification file.

        Returns
        -------
        Tuple[str, int]
            The PAM sequence and its signed size.

        Raises
        ------
        CrispritzPamError
            If the file cannot be read or declares a zero size.
        """
        # expected PAM file format:
        # NNNNNNNNNNNNNNNNNNNNNGG 3 (negative if pam upstream wrt grna)
        try:
            with open(fname, mode="r") as fin:
                pamseq, size = fin.readline().strip().split()
        except (IOError, Exception) as e:
            exception_handler(
                CrispritzPamError,
                f"Failed parsing PAM file {fname}",
                os.EX_IOERR,
                self._debug,
                e,
            )
        if size == 0:  # invalid pam size
            exception_handler(
                CrispritzPamError,
                f"Invalid PAM size: {size}",
                os.EX_DATAERR,
                self._debug,
            )
        return pamseq, int(size)

    def _refine_pam(self, pam_seq: str) -> str:
        """Return just the PAM motif from the full specification string.

        Slices the leading or trailing *size* bases of *pam_seq* according to
        the PAM orientation.

        Parameters
        ----------
        pam_seq : str
            The full PAM specification sequence.

        Returns
        -------
        str
            The isolated PAM motif.
        """
        # return just the pam sequence
        return pam_seq[: self._size] if self._upstream else pam_seq[-self._size :]

    def _assess_cas_system(self) -> None:
        """Classify the PAM motif against the known Cas-system catalogues.

        Sets the internal Cas-system identifier to the matching constant, or
        ``-1`` when the motif is not recognised.  Upstream and downstream
        motifs are matched against their respective catalogues.

        Returns
        -------
        None
        """
        self._cas_system = -1  # unknown cas system pam
        if self._upstream:
            if self._pam_seq in CPF1PAM:  # cpf1 cas system pam
                self._cas_system = CPF1
            elif self._pam_seq in CASXPAM:  # casx system pam
                self._cas_system = CASX
            elif self._pam_seq in SACAS9PAM:  # sacas9 system pam
                self._cas_system = SACAS9
        elif self._pam_seq in SPCAS9PAM:  # spcas9 system pam
            self._cas_system = SPCAS9
        elif self._pam_seq in XCAS9PAM:  # xcas9 pam
            self._cas_system = XCAS9

    @property
    def pamseq(self) -> str:
        """str: The isolated PAM motif."""
        return self._pam_seq

    @property
    def upstream(self) -> bool:
        """bool: ``True`` when the PAM lies upstream (5') of the spacer."""
        return self._upstream

    @property
    def size(self) -> int:
        """int: Length of the PAM motif."""
        return self._size

    @property
    def guide_size(self) -> int:
        """int: Spacer length implied by the specification (total minus PAM)."""
        return self._guide_size

    @property
    def cas_system(self) -> int:
        """int: Identifier of the recognised Cas system, or ``-1`` if unknown."""
        return self._cas_system
