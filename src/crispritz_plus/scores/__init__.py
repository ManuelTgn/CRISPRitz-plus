"""Off-target scoring subpackage for CRISPRitz-plus.

Re-exports :func:`~crispritz_plus.scores.shard_scoring.score_shards`, the entry
point that fills the ``cfd_score`` column of each per-partition target shard
using the CFD (Cutting Frequency Determination) model.  The CFD implementation
itself lives in the :mod:`~crispritz_plus.scores.cfd` subpackage.
"""

from .shard_scoring import score_shards
