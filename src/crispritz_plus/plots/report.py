"""
CLI driver for the ``generate-report`` subcommand in CRISPRitz-plus.

This module is a thin orchestration layer.  It reads and validates the
annotated search TSV, delegates all plot generation to
:mod:`.plots`, and emits structured progress output at the requested
verbosity level.

Pipeline overview
-----------------
1. Load and validate the annotated TSV via
   :func:`~.plots.load_annotated_targets`.
2. Detect annotation columns via
   :func:`~.plots.detect_annotation_columns` and log them.
3. Generate all figures via :func:`~.plots.generate_all_plots` and
   report the output paths.

Public API
----------
generate_report_cli
    Entry point called by the ``generate-report`` CLI subcommand.
"""

from time import time

import os


from ..crispritz_inputargs import CrispritzReportInputArgs
from ..exception_handlers import exception_handler
from ..verbosity import VERBOSITY_LVL, print_verbosity

from .crispritz_report_errors import CrispritzReportError
from .plots import detect_annotation_columns, load_annotated_targets, generate_all_plots


def generate_report_cli(args: CrispritzReportInputArgs) -> None:
    """Generate a graphical report from an annotated search output.

    This is the entry point for the ``generate-report`` CLI subcommand.
    It orchestrates the full report pipeline:

    1. **Load** - read and validate the annotated TSV via
       :func:`~.plots.load_annotated_targets`, converting the
       ``mismatches`` column to integers.
    2. **Detect** - identify annotation columns (those beyond the fixed
       search schema) and log them at verbose level.
    3. **Plot** - call :func:`~.plots.generate_all_plots` to produce one
       sequence-logo PNG per guide and one radar-chart PNG per
       ``(annotation_column, mismatch_level)`` combination per guide,
       writing all figures to ``args.outdir``.
    4. **Report** - log each output path at verbose level and emit a
       summary (figure count, output directory, elapsed time) at normal
       level.

    Parameters
    ----------
    args : CrispritzReportInputArgs
        Parsed CLI arguments.  The following attributes are consumed:

        * ``input_tsv`` - path to the annotated search-results TSV.
        * ``outdir`` - directory in which figures will be written.
        * ``mm`` - optional mismatch-level ceiling for the radar charts;
          when ``None`` the maximum observed level is used.
        * ``guides`` - optional explicit list of gap-stripped guide
          sequences to plot; when ``None`` all guides in the TSV are
          processed.
        * ``prefix`` - filename prefix applied to every output figure.
        * ``verbosity`` - integer controlling progress output volume
          (``0`` = silent, higher = more detail).
        * ``debug`` - when ``True``, exceptions propagate with a full
          traceback instead of a formatted user-facing message.

    Returns
    -------
    None
        Side-effect only: PNG figures are written to ``args.outdir``.

    Raises
    ------
    CrispritzReportError
        On TSV validation failure, empty input, missing guides, or any
        other error raised by the plot-generation pipeline.
    """
    start = time()
    print_verbosity(
        f"Generating graphical report from {args.input_tsv}",
        args.verbosity,
        VERBOSITY_LVL[1],  # Normal
    )
    try:
        df = load_annotated_targets(args.input_tsv, args.debug)
        annotation_columns = detect_annotation_columns(list(df.columns))
        print_verbosity(
            f"Detected {len(annotation_columns)} annotation column(s): "
            f"{', '.join(annotation_columns)}",
            args.verbosity,
            VERBOSITY_LVL[2],  # Verbose
        )
        outputs = generate_all_plots(
            tsv_path=args.input_tsv,
            outdir=args.outdir,
            debug=args.debug,
            mm=args.mm,
            guides=args.guides,
            prefix=args.prefix,
        )
    except CrispritzReportError as e:
        exception_handler(
            CrispritzReportError,
            "Report generation failed",
            os.EX_DATAERR,
            args.debug,
            e,
        )
    for path in outputs:
        print_verbosity(f"Wrote {path}", args.verbosity, VERBOSITY_LVL[2])
    print_verbosity(
        f"Report complete: {len(outputs)} figure(s) in {args.outdir} "
        f"({time() - start:.2f}s)",
        args.verbosity,
        VERBOSITY_LVL[1],
    )
