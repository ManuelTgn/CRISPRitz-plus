"""gRNA parsing and validation for CRISPRitz-plus.

Defines :class:`Guide`, a single validated guide RNA (with its reverse
complement), :class:`GuideList`, a length-consistent collection read from a
file, and the private :func:`_validate_guide_sequence` helper.  Guides are
supplied with the PAM region spelled out as ``N`` placeholders; parsing strips
that region and validates its placement against the :class:`~crispritz_plus.pam.PAM`
geometry.
"""

from typing import List

import os


from .crispritz_errors import CrispritzGuideError
from .dna_alphabet import reverse_complement, DNA
from .exception_handlers import exception_handler
from .pam import PAM


# ==============================================================================
# Guide class definition
# ==============================================================================


class Guide:
    """A single guide RNA sequence and its reverse complement.

    On construction the PAM placeholder region is stripped from *sequence*
    (validated to be all ``N``) and the reverse complement of the remaining
    spacer is computed and cached.

    Parameters
    ----------
    sequence : str
        The full guide sequence including the PAM placeholder region.
    pam : PAM
        The PAM describing the placeholder length and whether it is upstream.
    debug : bool
        When *True*, validation errors propagate with a full traceback.
    """

    def __init__(self, sequence: str, pam: PAM, debug: bool) -> None:
        self._debug = debug  # store debug flag
        # set guide sequence and its reverse complement
        self._sequence = _validate_guide_sequence(
            sequence, len(pam), pam.upstream, self._debug
        )
        self._sequence_rc = reverse_complement(self._sequence)

    def __len__(self) -> int:
        return len(self._sequence)

    @property
    def sequence(self) -> str:
        """str: The spacer sequence, with the PAM placeholder removed."""
        return self._sequence

    @property
    def reverse(self) -> str:
        """str: The reverse complement of the spacer sequence."""
        return self._sequence_rc


class GuideList:
    """A collection of equal-length guides read from a file.

    Reads one guide per line and enforces that every guide has the same
    length; both conditions are validated on construction.

    Parameters
    ----------
    fname : str
        Path to the guides file (one guide per line).
    pam : PAM
        PAM geometry applied when parsing each guide.
    debug : bool
        When *True*, parsing errors propagate with a full traceback.
    """

    def __init__(self, fname: str, pam: PAM, debug: bool) -> None:
        self._debug = debug  # store debug flag
        # read guides in input file
        self._guides = self._read_guides_file(pam, fname)

    def _read_guides_file(self, pam: PAM, fname: str) -> List[Guide]:
        """Parse the guides file into a list of :class:`Guide` objects.

        Parameters
        ----------
        pam : PAM
            PAM geometry applied when parsing each guide.
        fname : str
            Path to the guides file.

        Returns
        -------
        List[Guide]
            The parsed guides, all of equal length.

        Raises
        ------
        CrispritzGuideError
            If the file cannot be read, contains no guides, or holds guides of
            differing lengths.
        """
        # expected guides file format (2 guides):
        # CTAACAGTTGCTTTTATCACNNN
        # CTAACAGTTGCTTTTATCACNNN
        try:
            with open(fname, mode="r") as fin:
                guides = [Guide(line.strip(), pam, self._debug) for line in fin]
        except (IOError, Exception) as e:
            exception_handler(
                CrispritzGuideError,
                f"Failed parsing guide file {fname}",
                os.EX_IOERR,
                self._debug,
                e,
            )
        if not guides:  # no guide found in input file?
            exception_handler(
                CrispritzGuideError,
                f"No guide found in guides file {fname}",
                os.EX_IOERR,
                self._debug,
            )
        # guides must have same length
        if any(len(g) != len(guides[0]) for g in guides[1:]):
            exception_handler(
                CrispritzGuideError,
                "Mismatching guide lengths",
                os.EX_DATAERR,
                self._debug,
            )
        return guides

    @property
    def guides(self) -> List[Guide]:
        """List[Guide]: The parsed guides."""
        return self._guides


# ==============================================================================
# Internal helpers
# ==============================================================================


def _validate_guide_sequence(
    guide: str, pam_size: int, pam_upstream: bool, debug: bool
) -> str:
    """Validate the PAM placeholder and return the stripped spacer.

    Checks that the PAM-length slice of *guide* (leading slice when the PAM is
    upstream, trailing slice otherwise) is composed entirely of ``N``, then
    returns the guide with that slice removed.

    Parameters
    ----------
    guide : str
        The full guide sequence including the PAM placeholder region.
    pam_size : int
        Length of the PAM placeholder region.
    pam_upstream : bool
        ``True`` when the PAM lies upstream (5') of the spacer.
    debug : bool
        When *True*, validation errors propagate with a full traceback.

    Returns
    -------
    str
        The spacer sequence with the PAM placeholder removed.

    Raises
    ------
    CrispritzGuideError
        If the placeholder slice contains a non-``N`` base.
    """
    expected_pam_slice = guide[:pam_size] if pam_upstream else guide[-pam_size:]
    if any(nt != DNA[4] for nt in expected_pam_slice):
        exception_handler(
            CrispritzGuideError,
            f"Bad guide formatting, guide: {guide}",
            os.EX_DATAERR,
            debug,
        )
    return guide[pam_size:] if pam_upstream else guide[:-pam_size]
