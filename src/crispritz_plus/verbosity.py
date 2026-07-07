"""Verbosity-level constants and gated printing for CRISPRitz-plus.

Defines the four supported verbosity levels and a single helper,
:func:`print_verbosity`, that writes a message only when the caller's current
level meets a per-message threshold.  Every module in the package routes user
progress output through this helper so that the ``--verbosity`` CLI flag has a
uniform effect across the pipeline.

Verbosity levels
----------------
0
    Silent — suppress all progress output.
1
    Normal — high-level, once-per-stage messages.
2
    Detail — per-file / per-contig progress (single-threaded contexts).
3
    Debug — fine-grained tracing, safe to emit from inside parallel workers.

Module-level constants
----------------------
VERBOSITY_LVL : List[int]
    The ordered list of valid levels ``[0, 1, 2, 3]``.  Callers index into it
    (e.g. ``VERBOSITY_LVL[1]``) so that threshold references read symbolically
    rather than as bare integers.
"""

import sys


# ==============================================================================
# verbosity constant variables
# ==============================================================================

# enabled verbosity levels:
# 0 = silent
# 1 = normal
# 2 = detail
# 3 = debug
VERBOSITY_LVL = [0, 1, 2, 3]


# ==============================================================================
# verbosity functions
# ==============================================================================


def print_verbosity(message: str, verbosity: int, verbosity_threshold: int) -> None:
    """Print a message if the verbosity level meets the threshold.

    Writes the message to standard output if the current verbosity is greater
    than or equal to the specified threshold.

    Parameters
    ----------
    message:
        The message to print.
    verbosity:
        The current verbosity level.
    verbosity_threshold:
        The minimum verbosity level required to print the message.
    """
    if verbosity >= verbosity_threshold:
        sys.stdout.write(f"{message}\n")
    return
