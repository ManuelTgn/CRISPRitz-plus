"""CFD model loading and single-pair scoring (Doench et al., 2016).

Provides :func:`load_mismatch_pam_scores`, which unpickles the mismatch and PAM
penalty tables shipped alongside this module, and :func:`compute_cfd`, which
scores one guide / off-target pair against those tables.

Module-level constants
----------------------
MMSCORES : str
    Filename of the pickled mismatch score table (``"mismatch_score.pkl"``).
PAMSCORES : str
    Filename of the pickled PAM score table (``"pam_scores.pkl"``).
"""

from typing import Dict, Tuple

import pickle
import os


from crispritz_plus.dna_alphabet import dna2rna, reverse_complement
from crispritz_plus.exception_handlers import exception_handler

from ..crispritz_scores_errors import CfdScoreError


# ==============================================================================
# CFD score models
# ==============================================================================

MMSCORES = "mismatch_score.pkl"
PAMSCORES = "pam_scores.pkl"


# ==============================================================================
# Public API
# ==============================================================================


def load_mismatch_pam_scores(debug: bool) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Load the CFD mismatch and PAM penalty tables from disk.

    Unpickles the two model files shipped in the ``models`` directory next to
    this module (Doench et al., 2016).

    Parameters
    ----------
    debug : bool
        When *True*, load errors propagate with a full traceback instead of a
        formatted user-facing message.

    Returns
    -------
    Tuple[Dict[str, float], Dict[str, float]]
        The mismatch score table and the PAM score table, in that order.

    Raises
    ------
    CfdScoreError
        If either model file cannot be read.
    """
    modelspath = os.path.join(os.path.abspath(os.path.dirname(__file__)), "models")
    try:  # load mismatches and PAM scores (Doench et al., 2016)
        mmscores = pickle.load(open(os.path.join(modelspath, MMSCORES), mode="rb"))
        pamscores = pickle.load(open(os.path.join(modelspath, PAMSCORES), mode="rb"))
    except OSError as e:
        exception_handler(
            CfdScoreError,
            "An error occurred while loading CFD model files",
            os.EX_NOINPUT,
            debug,
            e,
        )
    return mmscores, pamscores


def compute_cfd(
    wildtype: str,
    sg: str,
    pam: str,
    mmscores: Dict[str, float],
    pamscores: Dict[str, float],
) -> float:  # sourcery skip: use-contextlib-suppress
    """Cutting Frequency Determination (CFD) score for one guide/off-target pair.

    The score is the product of per-position mismatch penalties (Doench et al.,
    2016) and a PAM-dinucleotide penalty, so a perfect on-target match scores
    ``1.0 * pamscores[pam]``.  Matched positions and bulge gaps (``'-'``)
    contribute a neutral factor of ``1.0``; the mismatch model is defined for
    substitutions only.

    Parameters
    ----------
    wildtype : str
        Guide body (PAM removed), 5'->3'.  Case-insensitive; ``T`` reads as
        ``U``.
    sg : str
        Off-target body (PAM removed), aligned to *wildtype* position by
        position.  Must be the same length as *wildtype*.
    pam : str
        PAM dinucleotide scored via *pamscores* (e.g. ``"GG"``).
    mmscores : Dict[str, float]
        Doench 2016 mismatch penalty table.
    pamscores : Dict[str, float]
        Doench 2016 PAM penalty table.

    Returns
    -------
    float
        CFD score in ``[0, 1]``.

    Raises
    ------
    ValueError
        If *wildtype* and *sg* differ in length (positional alignment is
        undefined).
    """
    if len(wildtype) != len(sg):
        raise ValueError(
            f"CFD requires equal-length guide and target; "
            f"got {len(wildtype)} and {len(sg)}"
        )
    # Normalise to uppercase RNA up front so the comparison is case-insensitive
    # (the target marks mismatches in lowercase) and 'T' is read as 'U'.
    wildtype = dna2rna(wildtype.upper())
    sg = dna2rna(sg.upper())
    score = 1.0
    for i, (wt_nt, ot_nt) in enumerate(zip(wildtype, sg)):
        if wt_nt == ot_nt:
            continue  # exact match: neutral factor
        if wt_nt == "-" or ot_nt == "-":
            if i == 0:
                return 0
            continue  # bulge gap: not penalised by the substitution model
        key = f"r{wt_nt}:d{reverse_complement(ot_nt)},{i + 1}"
        try:
            score *= mmscores[key]
        except KeyError:
            # No model entry (e.g. an ambiguous 'N' base, or a position beyond
            # the model's range): treat as neutral rather than abort the row.
            continue
    try:
        score *= pamscores[pam.upper()]
    except KeyError:
        # Non-canonical PAM dinucleotide: keep the mismatch-only product.
        pass
    return score
