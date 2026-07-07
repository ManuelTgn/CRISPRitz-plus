"""Off-target search subpackage for CRISPRitz-plus (``search``).

Exposes the CLI entry point :func:`search_offtargets_cli`, which unpacks the
validated search arguments and drives
:func:`~crispritz_plus.search.tst_explorer.search_offtargets_tst` — the
orchestrator that runs the per-partition C++ search, optional CFD scoring, and
the final sort/merge of the off-target table.
"""

from typing import List

from ..crispritz_inputargs import CrispritzSearchInputArgs
from .tst_explorer import search_offtargets_tst


def search_offtargets_cli(args: CrispritzSearchInputArgs) -> None:
    """Run the ``search`` subcommand.

    Unpacks the validated search arguments and delegates to
    :func:`~crispritz_plus.search.tst_explorer.search_offtargets_tst`.

    Parameters
    ----------
    args : CrispritzSearchInputArgs
        Validated CLI arguments (index partitions, PAM file, guides file,
        mismatch / bulge limits, output directory, threads, verbosity, debug
        flag, bulge mode, output mode, and the scoring flag).

    Returns
    -------
    None
    """
    search_offtargets_tst(
        args.indexes,
        args.pam_file,
        args.guides_file,
        args.mm,
        args.bdna,
        args.brna,
        args.outdir,
        args.threads,
        args.verbosity,
        args.debug,
        args.bulge_mode,
        args.output_mode,
        score=args.score,
    )


def search_offtargets(
    indexes: List[str],
    pam_file: str,
    guide_file: str,
    mm: int,
    bdna: int,
    brna: int,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
    bulge_mode: str = "mixed",
    output_mode: str = "both",
    sort_mode: str = "edit_distance",
    score: bool = False,
) -> None:
    """Exposed function to search off-targets on TST genome indexes.

    Delegates the off-targets search to
    :func:`~crispritz_plus.search.tst_explorer.search_offtargets_tst`.

    Parameters
    ----------
    indexes : List[str]
        Paths to the ``.bin`` TST partitions to search.
    pam_file : str
        Path to the PAM specification file.
    guides_file : str
        Path to the guides file (one guide per line).
    mm : int
        Maximum number of mismatches permitted in a hit.
    bdna : int
        Maximum number of DNA bulges permitted in a hit.
    brna : int
        Maximum number of RNA bulges permitted in a hit.
    outdir : str
        Directory for the output table, profiles, and temporary shards.
    threads : int
        Maximum number of worker threads.
    verbosity : int
        Verbosity level (see
        :data:`~crispritz_plus.verbosity.VERBOSITY_LVL`).
    debug : bool
        When *True*, errors propagate with a full traceback.
    bulge_mode : str, optional
        Whether DNA and RNA bulges may be combined. Defaults to ``"mixed"``.
    output_mode : str, optional
        Which result files to produce. Defaults to ``"both"``.
    sort_mode : str, optional
        Ordering applied to the merged table. Defaults to
        ``"edit_distance"``.
    score : bool, optional
        When ``True`` (and the PAM is SpCas9/xCas9), fill the CFD score column
        before merging. Defaults to ``False``.

    Returns
    -------
    None
    """
    search_offtargets_tst(
        indexes,
        pam_file,
        guide_file,
        mm,
        bdna,
        brna,
        outdir,
        threads,
        verbosity,
        debug,
        bulge_mode=bulge_mode,
        output_mode=output_mode,
        sort_mode=sort_mode,
        score=score,
    )
