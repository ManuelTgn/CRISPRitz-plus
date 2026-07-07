"""Validated CLI argument wrappers for every CRISPRitz-plus subcommand.

Each subcommand's parsed :class:`argparse.Namespace` is wrapped in a class that
validates paths, folders, and numeric ranges up front (via the shared
``_check_*`` / ``_validate_*`` helpers) and exposes the results as read-only
properties.  Validation failures are reported through the parser's ``error``
method, giving consistent usage-error messaging across the CLI.

Module-level constants
----------------------
OUTPUT_MODES : Tuple[str, ...]
    Canonical ``--output-mode`` choices, kept in sync with the C++
    ``OutputMode`` / ``OutputFormat`` enums and passed verbatim to the
    pybind11 layer.
BULGE_MODES : Tuple[str, ...]
    Canonical ``--bulge-mode`` choices, kept in sync with the C++
    ``BulgeMode`` enum.
"""

from argparse import Namespace
from glob import glob
from typing import List, Optional

import multiprocessing
import os


from .crispritz_argparse import CrispritzArgumentParser
from .verbosity import VERBOSITY_LVL


# Canonical CLI choices — kept in sync with the C++ SearchConfiguration enums
# (OutputMode / OutputFormat). These strings are passed verbatim to the
# pybind11 layer, which feeds them to output_mode_from_string() /
# output_format_from_string(); they must match those expected tokens exactly.
OUTPUT_MODES = ("both", "targets", "profile")

# Canonical bulge-mode CLI choices — kept in sync with the C++ BulgeMode enum.
# Passed verbatim to bulge_mode_from_string(): "mixed" allows DNA+RNA bulges in
# one alignment, "single" allows only one kind. (Legacy "both" still parses
# C++-side but is not offered as a CLI choice.)
BULGE_MODES = ("mixed", "single")


# ==============================================================================
# CLI argument wrapper classes
# ==============================================================================


class CrispritzInputArgs:
    """Base wrapper for validated subcommand arguments.

    Holds the parsed namespace and its parser, and provides the shared
    validators for the output folder, thread count, and verbosity that every
    subcommand reuses.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.
    parser : CrispritzArgumentParser
        The parser, used to report validation errors.
    """

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        self._args = args
        self._parser = parser

    def _validate_output_folder(self) -> None:
        """Validate the parent directory and create the output folder.

        Raises a usage error if the parent directory does not exist, then
        creates the output folder and stores its absolute path.

        Returns
        -------
        None
        """
        outdir = os.path.abspath(self._args.outdir)
        parentdir = os.path.dirname(self._args.outdir) or os.getcwd()
        _check_folder(
            parentdir, self._parser, f"Cannot create output folder {self._args.outdir}"
        )
        os.makedirs(outdir, exist_ok=True)  # create output folder
        self._outdir = outdir

    def _validate_threads(self) -> None:
        """Validate and store the requested thread count.

        Returns
        -------
        None
        """
        self._threads = _validate_threads_num(self._args.threads, self._parser)

    def _validate_verbosity(self) -> None:
        """Validate and store the requested verbosity level.

        Returns
        -------
        None
        """
        self._verbosity = _validate_verbosity_value(self._args.verbosity, self._parser)

    @property
    def outdir(self) -> str:
        """str: Absolute path of the validated output folder."""
        return self._outdir

    @property
    def threads(self) -> int:
        """int: Requested number of threads."""
        return self._args.threads

    @property
    def verbosity(self) -> int:
        """int: Requested verbosity level."""
        return self._args.verbosity

    @property
    def debug(self) -> bool:
        """bool: Whether debug mode is enabled."""
        return self._args.debug


class CrispritzEnrichmentInputArgs(CrispritzInputArgs):
    """Validated arguments for the ``add-variants`` subcommand.

    Extends :class:`CrispritzInputArgs` with VCF- and genome-folder discovery
    and the indel / keep flags, validating all inputs on construction.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.
    parser : CrispritzArgumentParser
        The parser used to report validation errors.
    """

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()

    def _validate_vcf_folder(self) -> None:
        """Validate the VCF folder and collect its ``*.vcf.gz`` files.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            Via the parser's ``error`` method if the folder is missing or
            contains no VCF files.
        """
        _check_folder(
            self._args.vcf, self._parser, f"Cannot find VCF folder {self._args.vcf}"
        )
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        _check_retrieved_files(
            self._vcfs, self._parser, f"No VCF file found in {self._args.vcf}"
        )

    def _validate_genome_folder(self) -> None:
        """Validate the genome folder and collect its FASTA files.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            Via the parser's ``error`` method if the folder is missing or
            contains no ``*.fa`` / ``*.fasta`` files.
        """
        _check_folder(
            self._args.genome,
            self._parser,
            f"Cannot find input genome folder {self._args.genome}",
        )
        self._fastas = glob(os.path.join(self._args.genome, "*.fa")) + glob(
            os.path.join(self._args.genome, "*.fasta")
        )
        _check_retrieved_files(
            self._fastas, self._parser, f"No FASTA file found in {self._args.genome}"
        )

    def _check_consistency(self) -> None:
        """Run all enrichment input validations in order.

        Returns
        -------
        None
        """
        self._validate_vcf_folder()  # check vcf folder
        self._validate_genome_folder()  # check genome folder
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    @property
    def vcfs(self) -> List[str]:
        """List[str]: Discovered VCF file paths."""
        return self._vcfs

    @property
    def fastas(self) -> List[str]:
        """List[str]: Discovered FASTA file paths."""
        return self._fastas

    @property
    def indels(self) -> bool:
        """bool: Whether indel processing is enabled."""
        return self._args.indels

    @property
    def keep(self) -> bool:
        """bool: Whether all variants are kept regardless of FILTER."""
        return self._args.keep


class CrispritzIndexingInputArgs(CrispritzInputArgs):
    """Validated arguments for the ``index-genome`` subcommand.

    Extends :class:`CrispritzInputArgs` with genome-folder discovery, PAM-file
    and maximum-bulge validation.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.
    parser : CrispritzArgumentParser
        The parser used to report validation errors.
    """

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()  # check input args consistency

    def _validate_genome_folder(self) -> None:
        """Validate the genome folder and collect its FASTA files.

        Returns
        -------
        None
        """
        _check_folder(
            self._args.genome,
            self._parser,
            f"Cannot find input genome folder {self._args.genome}",
        )
        self._fastas = glob(os.path.join(self._args.genome, "*.fa")) + glob(
            os.path.join(self._args.genome, "*.fasta")
        )
        _check_retrieved_files(
            self._fastas, self._parser, f"No FASTA file found in {self._args.genome}"
        )

    def _validate_pam_file(self) -> None:
        """Validate that the PAM file exists.

        Returns
        -------
        None
        """
        _check_file(
            self._args.pam_file,
            self._parser,
            f"Cannot find input PAM file {self._args.pam_file}",
        )

    def _validate_bmax(self) -> None:
        """Validate that the maximum-bulge value is non-negative.

        Returns
        -------
        None
        """
        if self._args.bmax < 0:
            self._parser.error(f"Invalid max bulge value: {self._args.bmax}")

    def _check_consistency(self) -> None:
        """Run all indexing input validations in order.

        Returns
        -------
        None
        """
        self._validate_genome_folder()  # check genome folder
        self._validate_pam_file()  # check pam file
        self._validate_bmax()  # check max bulge
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    @property
    def fastas(self) -> List[str]:
        """List[str]: Discovered FASTA file paths."""
        return self._fastas

    @property
    def pam_file(self) -> str:
        """str: Path to the PAM specification file."""
        return self._args.pam_file

    @property
    def bmax(self) -> int:
        """int: Maximum number of bulges to index."""
        return self._args.bmax


class CrispritzSearchInputArgs(CrispritzInputArgs):
    """Validated arguments for the ``search`` subcommand.

    Extends :class:`CrispritzInputArgs` with index-folder discovery, PAM and
    guides file checks, and mismatch / bulge range validation.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.
    parser : CrispritzArgumentParser
        The parser used to report validation errors.
    """

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()

    def _check_consistency(self) -> None:
        """Run all search input validations in order.

        Returns
        -------
        None
        """
        self._validate_index_genome()  # check genome index folder
        self._validate_pam_file()  # check pam file
        self._validate_guides_file()  # check guides file
        self._validate_mm()  # check mismatch
        self._validate_bdna()  # check dna bulge
        self._validate_brna()  # check rna bulge
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    def _validate_index_genome(self) -> None:
        """Validate the index folder and collect its ``*.bin`` partitions.

        Returns
        -------
        None
        """
        _check_folder(
            self._args.genome_index,
            self._parser,
            f"Cannot find genome index folder {self._args.genome_index}",
        )
        self._indexes = glob(os.path.join(self._args.genome_index, "*.bin"))
        _check_retrieved_files(
            self._indexes,
            self._parser,
            f"No TST index found in {self._args.genome_index}",
        )

    def _validate_pam_file(self) -> None:
        """Validate that the PAM file exists.

        Returns
        -------
        None
        """
        _check_file(
            self._args.pam_file,
            self._parser,
            f"Cannot find input PAM file {self._args.pam_file}",
        )

    def _validate_guides_file(self) -> None:
        """Validate that the guides file exists.

        Returns
        -------
        None
        """
        _check_file(
            self._args.guides_file,
            self._parser,
            f"Cannot find input guides file {self._args.guides_file}",
        )

    def _validate_mm(self) -> None:
        """Validate that the mismatch count is non-negative.

        Returns
        -------
        None
        """
        if self._args.mm < 0:
            self._parser.error(f"Invalid max mismatch value: {self._args.mm}")

    def _validate_bdna(self) -> None:
        """Validate that the maximum DNA-bulge count is non-negative.

        Returns
        -------
        None
        """
        if self._args.bdna < 0:
            self._parser.error(f"Invalid max DNA bulge value: {self._args.bdna}")

    def _validate_brna(self) -> None:
        """Validate that the maximum RNA-bulge count is non-negative.

        Returns
        -------
        None
        """
        if self._args.brna < 0:
            self._parser.error(f"Invalid max RNA bulge value: {self._args.brna}")

    @property
    def indexes(self) -> List[str]:
        """List[str]: Discovered ``*.bin`` TST partition paths."""
        return self._indexes

    @property
    def pam_file(self) -> str:
        """str: Path to the PAM specification file."""
        return self._args.pam_file

    @property
    def guides_file(self) -> str:
        """str: Path to the guides file."""
        return self._args.guides_file

    @property
    def mm(self) -> int:
        """int: Maximum number of mismatches allowed."""
        return self._args.mm

    @property
    def bdna(self) -> int:
        """int: Maximum number of DNA bulges allowed."""
        return self._args.bdna

    @property
    def brna(self) -> int:
        """int: Maximum number of RNA bulges allowed."""
        return self._args.brna

    @property
    def output_mode(self) -> str:
        """str: Selected output mode (see :data:`OUTPUT_MODES`)."""
        return self._args.output_mode

    @property
    def bulge_mode(self) -> str:
        """str: Selected bulge mode (see :data:`BULGE_MODES`)."""
        return self._args.bulge_mode

    @property
    def score(self) -> bool:
        """bool: Whether off-target scoring is enabled."""
        return self._args.score


class CrispritzAnnotateInputArgs(CrispritzInputArgs):
    """Validated arguments for the ``annotate-results`` subcommand.

    Extends :class:`CrispritzInputArgs` with targets-file, annotation-BED, and
    annotation-name validation.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.
    parser : CrispritzArgumentParser
        The parser used to report validation errors.
    """

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()

    def _check_consistency(self) -> None:
        """Run all annotation input validations in order.

        Returns
        -------
        None
        """
        self._validate_targets_file()  # check search targets table
        self._validate_annotations()  # check annotation BED files
        self._validate_annotation_names()  # check names/files length match
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    def _validate_targets_file(self) -> None:
        """Validate that the search targets file exists.

        Returns
        -------
        None
        """
        _check_file(
            self._args.targets_file,
            self._parser,
            f"Cannot find input targets file {self._args.targets_file}",
        )

    def _validate_annotations(self) -> None:
        """Validate the annotation BED files exist and have a supported suffix.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            Via the parser's ``error`` method if no annotation is provided, a
            file is missing, or a file is not ``.bed`` / ``.bed.gz``.
        """
        if not self._args.annotations:
            self._parser.error("No annotation BED file provided")
        for bed in self._args.annotations:
            _check_file(bed, self._parser, f"Cannot find annotation file {bed}")
            if not bed.endswith((".bed", ".bed.gz")):
                self._parser.error(
                    f"Unsupported annotation file '{bed}'. "
                    "Expected a .bed or .bed.gz file"
                )

    def _validate_annotation_names(self) -> None:
        """Validate that annotation names match the number of annotation files.

        Returns
        -------
        None
        """
        names = self._args.annotation_names
        if names is not None and len(names) != len(self._args.annotations):
            self._parser.error(
                f"Number of --annotation-names ({len(names)}) does not match "
                f"the number of --annotations files ({len(self._args.annotations)})"
            )

    @property
    def targets_file(self) -> str:
        """str: Path to the search targets TSV."""
        return self._args.targets_file

    @property
    def annotations(self) -> List[str]:
        """List[str]: Annotation BED file paths."""
        return self._args.annotations

    @property
    def annotation_names(self) -> Optional[List[str]]:
        """Optional[List[str]]: Custom annotation column names, or ``None``."""
        return self._args.annotation_names


class CrispritzReportInputArgs(CrispritzInputArgs):
    """Validated arguments for the generate-report subcommand.

    Extends :class:`CrispritzInputArgs` with input targets-file, number of
    mismatches, target guide, and report prefix name.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.
    parser : CrispritzArgumentParser
        The parser used to report validation errors.
    """

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._validate()

    def _validate(self) -> None:
        """Validate the input TSV, mismatch level, and verbosity.

        Returns
        -------
        None

        Raises
        ------
        SystemExit
            Via the parser's ``error`` method if the input TSV is missing or
            the mismatch level is negative.
        """
        _check_file(
            self._args.input_tsv,
            self._parser,
            f"Cannot find annotated results file {self._args.input_tsv}",
        )
        if self._args.mm is not None and self._args.mm < 0:
            self._parser.error(f"Invalid mismatch level: {self._args.mm}")

    @property
    def input_tsv(self) -> str:
        """str: Path to the annotated results TSV."""
        return self._args.input_tsv

    @property
    def mm(self):
        """Optional[int]: Maximum mismatch level to plot, or ``None``."""
        return self._args.mm

    @property
    def guides(self):
        """Optional[Sequence[str]]: Explicit guides to plot, or ``None``."""
        return self._args.guides

    @property
    def prefix(self) -> str:
        """str: Output filename prefix."""
        return self._args.prefix


# ==============================================================================
# Internal helpers
# ==============================================================================


def _check_folder(dirname: str, parser: CrispritzArgumentParser, msg: str) -> None:
    """Report a usage error if *dirname* is not an existing directory.

    Parameters
    ----------
    dirname : str
        Directory path to check.
    parser : CrispritzArgumentParser
        Parser used to report the error.
    msg : str
        Error message shown on failure.

    Returns
    -------
    None
    """
    if not os.path.exists(dirname) or not os.path.isdir(dirname):
        parser.error(msg)


def _check_file(fname: str, parser: CrispritzArgumentParser, msg: str) -> None:
    """Report a usage error if *fname* is not an existing file.

    Parameters
    ----------
    fname : str
        File path to check.
    parser : CrispritzArgumentParser
        Parser used to report the error.
    msg : str
        Error message shown on failure.

    Returns
    -------
    None
    """
    if not os.path.exists(fname) or not os.path.isfile(fname):
        parser.error(msg)


def _check_retrieved_files(
    fnames: List[str], parser: CrispritzArgumentParser, msg: str
) -> None:
    """Report a usage error if *fnames* is empty.

    Parameters
    ----------
    fnames : List[str]
        Collected file paths.
    parser : CrispritzArgumentParser
        Parser used to report the error.
    msg : str
        Error message shown on failure.

    Returns
    -------
    None
    """
    if not fnames:
        parser.error(msg)


def _validate_threads_num(threads: int, parser: CrispritzArgumentParser) -> int:
    """Validate a thread count against the available CPU cores.

    A value of ``0`` is interpreted as "use all cores".

    Parameters
    ----------
    threads : int
        Requested number of threads.
    parser : CrispritzArgumentParser
        Parser used to report the error.

    Returns
    -------
    int
        The validated thread count (all cores when *threads* is ``0``).

    Raises
    ------
    SystemExit
        Via the parser's ``error`` method if *threads* is negative or exceeds
        the available cores.
    """
    max_threads = multiprocessing.cpu_count()
    if threads < 0 or threads > max_threads:
        parser.error(
            f"Forbidden number of threads provided ({threads}). "
            f"Max number of available cores: {max_threads}"
        )
    return max_threads if threads == 0 else threads


def _validate_verbosity_value(verbosity: int, parser: CrispritzArgumentParser) -> int:
    """Validate a verbosity level against :data:`~crispritz_plus.verbosity.VERBOSITY_LVL`.

    Parameters
    ----------
    verbosity : int
        Requested verbosity level.
    parser : CrispritzArgumentParser
        Parser used to report the error.

    Returns
    -------
    int
        The validated verbosity level.

    Raises
    ------
    SystemExit
        Via the parser's ``error`` method if the level is not recognised.
    """
    if verbosity not in VERBOSITY_LVL:
        parser.error(f"Forbidden verbosity level selected ({verbosity})")
    return verbosity
