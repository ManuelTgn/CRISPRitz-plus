"""Shared filesystem and process helpers for CRISPRitz-plus.

Small, dependency-light utilities used across the package: directory / file
creation and validation, index-file presence checks, worker-count sizing, and
file rename / removal.  Also defines the tool naming and subcommand constants.

Module-level constants
----------------------
TOOLNAME : str
    Display name of the tool (``"CRISPRitz"``).
COMMAND : str
    The console command used to invoke the tool (``"crispritz"``).
SUBCOMMANDS : List[str]
    The ordered list of CLI subcommands exposed by the tool.
"""

from pathlib import Path

import os


# ==============================================================================
# Tool-wide constants
# ==============================================================================

#: Human-facing display name of the tool.
TOOLNAME = "CRISPRitz"

#: Console command used to invoke the tool.
COMMAND = "crispritz"

#: Ordered list of CLI subcommands exposed by the tool.
SUBCOMMANDS = [
    "add-variants",
    "index-genome",
    "search",
    "annotate-results",
    "generate-report",
]


# ==============================================================================
# Tool-wide utils functions
# ==============================================================================


def validate_directory(path: str, create: bool = False) -> Path:
    """Resolve and validate a directory path, optionally creating it.

    Parameters
    ----------
    path : str
        Directory path to validate.
    create : bool, optional
        When *True*, create the directory (and parents) if it does not exist.
        Defaults to ``False``.

    Returns
    -------
    Path
        The resolved absolute directory path.

    Raises
    ------
    FileNotFoundError
        If the directory does not exist and *create* is ``False``.
    NotADirectoryError
        If the resolved path exists but is not a directory.
    """
    dir_path = Path(path).resolve()
    if not dir_path.exists():
        if create:
            dir_path.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Directory not found: {path}")
    if not dir_path.is_dir():
        raise NotADirectoryError(f"Not a directory: {path}")
    return dir_path


def create_folder(folder: str) -> str:
    """Create a folder if needed and return its absolute path.

    Parameters
    ----------
    folder : str
        Path of the folder to create; existing folders are accepted.

    Returns
    -------
    str
        The absolute path of the folder.
    """
    os.makedirs(folder, exist_ok=True)
    assert os.path.isdir(folder)
    return os.path.abspath(folder)


def find_fasta_index(fasta_path: str) -> bool:
    """Return whether a FASTA ``.fai`` index exists for *fasta_path*.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file (without the ``.fai`` suffix).

    Returns
    -------
    bool
        ``True`` if ``{fasta_path}.fai`` exists, ``False`` otherwise.
    """
    return os.path.isfile(f"{fasta_path}.fai")


def find_tabix_index(fname_path: str) -> bool:
    """Return whether a tabix ``.tbi`` index exists for *fname_path*.

    Parameters
    ----------
    fname_path : str
        Path to the compressed file (without the ``.tbi`` suffix).

    Returns
    -------
    bool
        ``True`` if ``{fname_path}.tbi`` exists, ``False`` otherwise.
    """
    return os.path.isfile(f"{fname_path}.tbi")


def set_processes(tasks_n: int, threads: int) -> int:
    """Return the number of worker processes to use.

    Caps the requested thread count at the number of available tasks so that
    no idle workers are spawned.

    Parameters
    ----------
    tasks_n : int
        Number of independent tasks to process.
    threads : int
        Maximum number of worker processes requested.

    Returns
    -------
    int
        ``min(tasks_n, threads)``.
    """
    return min(tasks_n, threads)


def rename_files(origin: str, dest: str) -> str:
    """Rename / move a file and return the destination path.

    Parameters
    ----------
    origin : str
        Existing source path.
    dest : str
        Destination path.

    Returns
    -------
    str
        The destination path *dest*.
    """
    origin_path = Path(origin)
    origin_path.rename(dest)
    assert os.path.isfile(dest)
    return dest


def remove_file(path: str) -> None:
    """Delete a file from disk.

    Parameters
    ----------
    path : str
        Path of the file to remove.

    Returns
    -------
    None
    """
    Path(path).unlink()
