from ..crispritz_inputargs import CrispritzAnnotateInputArgs
from .annotation import annotate_results


def annotate_results_cli(args: CrispritzAnnotateInputArgs) -> None:
    annotate_results(
        args.targets_file,
        args.annotations,
        args.outdir,
        args.annotation_names,
        args.threads,
        args.verbosity,
        args.debug,
    )