"""Custom ``argparse`` parser and help formatting for CRISPRitz-plus.

Defines :class:`CrispritzArgumentParser`, an :class:`argparse.ArgumentParser`
subclass that renders errors in red, substitutes the tool version into usage
strings, and installs a help formatter that honours ``SUPPRESS`` on usage.
Centralising these tweaks keeps the CLI's look and error behaviour consistent
across every subcommand.
"""

from argparse import (
    SUPPRESS,
    ArgumentParser,
    HelpFormatter,
    Action,
    _MutuallyExclusiveGroup,
)
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn
from colorama import Fore

import sys
import os


from .utils import COMMAND
from .version import __version__


# define abstract generic types for typing
_D = TypeVar("_D")
_V = TypeVar("_V")


class CrispritzArgumentParser(ArgumentParser):
    """An ``ArgumentParser`` with CRISPRitz-plus help and error styling.

    Customises three behaviours relative to :class:`argparse.ArgumentParser`:
    a nested help formatter that suppresses usage when requested, substitution
    of the tool version into ``usage`` strings, and coloured error reporting
    that points the user at ``-h``.
    """

    class CrispritzHelpFormatter(HelpFormatter):
        """Help formatter that honours ``SUPPRESS`` for the usage section."""

        def add_usage(  # type: ignore
            self,
            usage: str,
            actions: Iterable[Action],
            groups: Iterable[_MutuallyExclusiveGroup],
            prefix: Optional[str] = None,
        ) -> None:
            """Add the usage section unless it is suppressed.

            Parameters
            ----------
            usage : str
                The usage string, or :data:`argparse.SUPPRESS` to omit it.
            actions : Iterable[Action]
                The parser actions to include in the usage line.
            groups : Iterable[_MutuallyExclusiveGroup]
                Mutually-exclusive argument groups.
            prefix : Optional[str], optional
                Usage prefix. Defaults to ``None``.

            Returns
            -------
            None
            """
            if usage != SUPPRESS:
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # initialize the formatter

    def __init__(self, *args: Tuple[_D], **kwargs: Dict[_D, _V]) -> None:
        # set custom help formatter defined as
        kwargs["formatter_class"] = self.CrispritzHelpFormatter  # type: ignore
        # replace the default version display in usage help with a custom
        # version display formatter
        if "usage" in kwargs:
            kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)  # type: ignore
        # initialize argument parser object with input parameters for
        # usage display
        super().__init__(*args, **kwargs)  # type: ignore

    def error(self, error: str) -> NoReturn:  # type: ignore
        """Report a usage error in red and exit.

        Writes a coloured ``ERROR:`` message and a pointer to ``-h`` to
        standard error, then exits with :data:`os.EX_USAGE`.

        Parameters
        ----------
        error : str
            The error message to display.

        Returns
        -------
        NoReturn
            Does not return; the process exits.
        """
        errormsg = (
            f"{Fore.RED}\nERROR: {error}.{Fore.RESET}"
            + f"\n\nRun {COMMAND} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error

    def error_noargs(self) -> None:
        """Print full help and exit when invoked with no arguments.

        Returns
        -------
        None
            Does not return; the process exits with :data:`os.EX_NOINPUT`.
        """
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code
