"""Ternary Search Tree index construction for CRISPRitz-plus.

Provides :func:`build_ternary_search_tree`, the Python orchestration layer for
the ``index-genome`` stage.  It reads each per-chromosome FASTA, derives the
contig name, and hands the sequence and PAM geometry to the C++ builder, which
writes the nibble-packed ``.bin`` partition files that the search stage loads.
"""

from typing import List

import os


from ..crispritz_cpp import build_tree_cpp
from ..crispritz_errors import CrispritzTstError
from ..exception_handlers import exception_handler
from ..genome_io import GenomeReader
from ..pam import PAM
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

    Reads each FASTA, extracts the chromosome name from the header, and calls
    the C++ TST builder.  One or more ``.bin`` partition files are written per
    chromosome into *outdir*.

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

    Raises
    ------
    CrispritzTstError
        If the C++ builder fails on any input FASTA.
    """
    pam = PAM(pam_file, debug)
    for fasta in fastas:
        reader = GenomeReader(fasta, debug)
        reader.read()
        # the contig name is used in the output .bin filename(s).
        # Add leading 'chr' prefix to improve the legacy naming used by
        # (even though the search binary expects e.g. "1" not "chr1").
        contig = (
            reader.header if reader.header.startswith("chr") else f"chr{reader.header}"
        )
        print_verbosity(
            f"Building TST index for {contig} ({fasta})", verbosity, VERBOSITY_LVL[2]
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
            )
        except Exception as e:
            exception_handler(
                CrispritzTstError,
                f"Failed building ternary search tree on {fasta}",
                os.EX_DATAERR,
                debug,
                e,
            )
