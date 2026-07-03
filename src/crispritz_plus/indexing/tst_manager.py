""" """

from time import time
from typing import List

import os


from ..crispritz_cpp import build_tree_cpp
from ..crispritz_errors import CrispritzTstError
from ..exception_handlers import exception_handler
from ..genome_io import GenomeReader
from ..pam import PAM
from ..progress import progress_bar
from ..verbosity import VERBOSITY_LVL, print_verbosity


def build_ternary_search_tree(
    fastas: List[str],
    pam_file: str,
    bmax: int,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """Build a Ternary Search Tree index for every input FASTA file.

    Reads each FASTA, extracts the chromosome name from the header, and
    calls the C++ TST builder.  One or more ``.bin`` partition files are
    written per chromosome into the current working directory.

    Parameters
    ----------
    fastas:
        Paths to the per-chromosome FASTA files to index.
    pam_file:
        Path to the PAM specification file.
    bmax:
        Maximum number of bulges; passed to the C++ builder so it extracts
        enough extra bases per site to support bulge-aware search.
    outdir:
        Path to the directory where the genome index will be stored.
    threads:
        Number of OpenMP threads for the PAM search phase.
    verbosity:
        Verbosity level (see ``VERBOSITY_LVL``).
    debug:
        When *True*, exceptions propagate with full stack traces.
    """
    pam = PAM(pam_file, debug)
    print_verbosity(
        f"build_ternary_search_tree: pam={pam.pamseq}, pam_size={pam.size}, "
        f"guide_size={pam.guide_size}, upstream={pam.upstream}, bmax={bmax}, "
        f"threads={threads}, outdir={outdir!r}",
        verbosity,
        VERBOSITY_LVL[3],
    )
    print_verbosity(
        f"Building genome index for {len(fastas)} FASTA file(s)",
        verbosity,
        VERBOSITY_LVL[1],
    )
    start = time()  # track total time
    for i, fasta in enumerate(progress_bar(fastas, "Constructed TST indexes", verbosity), start=1):
        reader = GenomeReader(fasta, debug)
        reader.read()
        # the contig name is used in the output .bin filename(s).
        # Add leading 'chr' prefix to improve the legacy naming used by
        # (even though the search binary expects e.g. "1" not "chr1").
        contig = (
            reader.header if reader.header.startswith("chr") else f"chr{reader.header}"
        )
        print_verbosity(
            f"[{i}/{len(fastas)}] Building TST index for {contig} ({fasta})",
            verbosity,
            VERBOSITY_LVL[2],
        )
        try:
            build_tree_cpp(
                reader.to_string(),
                contig,
                pam.pamseq,
                pam.guide_size + pam.size,
                pam.size,
                pam.upstream,
                outdir,
                bmax,
                threads,
                verbosity,
            )
        except Exception as e:
            exception_handler(
                CrispritzTstError,
                f"Failed building ternary search tree on {fasta}",
                os.EX_DATAERR,
                debug,
                e,
            )
    print_verbosity(
        f"Genome index built for {len(fastas)} FASTA file(s) in "
        f"{time() - start:.2f}s",
        verbosity,
        VERBOSITY_LVL[1],
    )
