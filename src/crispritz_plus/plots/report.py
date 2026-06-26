""" """

from time import time

import os


from ..crispritz_inputargs import CrispritzReportInputArgs
from ..exception_handlers import exception_handler
from ..verbosity import VERBOSITY_LVL, print_verbosity

from .crispritz_report_errors import CrispritzReportError
from .plots import detect_annotation_columns, load_annotated_targets, generate_all_plots


def generate_report_cli(args: CrispritzReportInputArgs) -> None:
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
