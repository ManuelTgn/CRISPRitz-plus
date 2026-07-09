"""Centralised error handling and signal handling for CRISPRitz-plus.

Provides the two entry points every module uses to fail consistently:
:func:`exception_handler`, which either re-raises with a full traceback (debug
mode) or prints a coloured error and exits with a specific code, and
:func:`sigint_handler`, which exits gracefully on keyboard interrupt.

Routing all failures through :func:`exception_handler` keeps error formatting,
exit codes, and the debug/production behaviour split identical across the
whole package.
"""

from typing import Optional, NoReturn
from colorama import init, Fore

import sys
import os


def sigint_handler() -> None:
    """Handle SIGINT (keyboard interrupt) by exiting gracefully.

    Writes a short notice to standard error and terminates the process with
    :data:`os.EX_OSERR`.

    Returns
    -------
    None
        Does not return; the process exits.
    """
    # print message when SIGINT is caught to exit gracefully from the execution
    sys.stderr.write(f"\nCaught SIGINT. Exit CRISPRitz+\n")
    sys.exit(os.EX_OSERR)  # mark as os error code


def exception_handler(
    exception_type: type,
    exception: str,
    code: int,
    debug: bool,
    e: Optional[Exception] = None,
) -> NoReturn:
    """Raise or report an error, then exit, honouring debug mode.

    In debug mode the error is raised as *exception_type* with a full
    traceback (chained from *e* when supplied), so the stack is preserved for
    diagnosis.  Otherwise a red ``ERROR:`` line is written to standard error
    and the process exits with *code*.

    Parameters
    ----------
    exception_type : type
        The exception class to raise in debug mode.
    exception : str
        The human-readable error message.
    code : int
        Process exit code used in non-debug mode (typically an ``os.EX_*``
        constant).
    debug : bool
        When *True*, raise with a full traceback instead of exiting.
    e : Optional[Exception], optional
        A previous exception to chain from (``raise ... from e``).  Defaults
        to ``None``.

    Returns
    -------
    NoReturn
        Never returns normally: it either raises (debug) or exits the process.
    """
    init()  # initialize colorama render
    if debug:  # debug mode -> always trace back the full error stack
        if e:  # inherits from previous error
            raise exception_type(f"\n\n{exception}") from e
        raise exception_type(f"\n\n{exception}")  # divide exception message from stack
    # gracefully trigger error and exit execution
    sys.stderr.write(f"{Fore.RED}\n\nERROR: {exception}\n{Fore.RESET}")
    sys.exit(code)  # exit execution returning appropriate error code
