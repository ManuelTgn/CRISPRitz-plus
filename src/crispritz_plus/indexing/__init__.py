"""Genome-indexing subpackage for CRISPRitz-plus (``index-genome``).

Exposes the CLI entry point :func:`index_genome_cli`, which unpacks the
validated indexing arguments and drives
:func:`~crispritz_plus.indexing.tst_manager.build_ternary_search_tree` to
construct the Ternary Search Tree (TST) index consumed by the ``search``
stage.
"""

from typing import List


from ..crispritz_inputargs import CrispritzIndexingInputArgs
from .tst_manager import build_ternary_search_tree


def index_genome_cli(args: CrispritzIndexingInputArgs) -> None:
    """Run the ``index-genome`` subcommand.

    Unpacks the validated indexing arguments and delegates to
    :func:`~crispritz_plus.indexing.tst_manager.build_ternary_search_tree`.

    Parameters
    ----------
    args : CrispritzIndexingInputArgs
        Validated CLI arguments. The FASTA list, PAM file, maximum-bulge
        value, output directory, thread count, verbosity, and debug flag are
        forwarded to the builder.

    Returns
    -------
    None
    """
    build_ternary_search_tree(
        args.fastas,
        args.pam_file,
        args.bmax,
        args.outdir,
        args.threads,
        args.verbosity,
        args.debug,
    )


def index_genome(
    fastas: List[str],
    pam_file: str,
    bmax: int,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """Exposed function to build a TST genome index.

    Delegates the TST construction to
    :func:`~crispritz_plus.indexing.tst_manager.build_ternary_search_tree`.

    Parameters
    ----------
    fastas : List[str]
        Paths to the per-chromosome FASTA files to index.
    pam_file : str
        Path to the PAM specification file.
    bmax : int
        Maximum number of bulges; passed to the C++ builder so it extracts
        enough extra bases per site to support bulge-aware search.
    outdir : str
        Directory where the genome index will be stored.
    threads : int
        Number of OpenMP threads for the PAM search phase.
    verbosity : int
        Verbosity level (see
        :data:`~crispritz_plus.verbosity.VERBOSITY_LVL`).
    debug : bool
        When *True*, exceptions propagate with full stack traces.

    Returns
    -------
    None
    """
    build_ternary_search_tree(fastas, pam_file, bmax, outdir, threads, verbosity, debug)
