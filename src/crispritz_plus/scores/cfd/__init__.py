"""CFD (Cutting Frequency Determination) scoring for CRISPRitz-plus.

Implements the Doench et al. (2016) CFD off-target model.  Re-exports
:func:`~crispritz_plus.scores.cfd.cfdscore.load_mismatch_pam_scores` (loads the
mismatch and PAM penalty tables) and
:func:`~crispritz_plus.scores.cfd.cfdscore.compute_cfd` (scores a single
guide / off-target pair).
"""

from .cfdscore import load_mismatch_pam_scores, compute_cfd
