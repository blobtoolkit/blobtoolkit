#!/usr/bin/env python3

"""
BlobToolKit - assembly exploration, QC and filtering.

usage: btk [<command>] [<subcommand>] [<args>...] [-h|--help] [--version]

commands:
    blobtools       BTK command line component
    pipeline        BTK analysis pipeline
    -h, --help      show this
    -v, --version   show version number
See 'btk <command> --help' for more information on a specific command.

"""


import os
import sys

from docopt import DocoptExit
from docopt import docopt
from pkg_resources import working_set
from tolkein import tolog

from .lib.version import __version__

LOGGER = tolog.logger(__name__)


def cli():
    """Entry point."""
    if len(sys.argv) > 2:
        try:
            args = docopt(__doc__, help=False, version=__version__)
        except DocoptExit:
            args = {
                "<command>": " ".join(sys.argv[1:2]),
            }
    elif len(sys.argv) > 1:
        try:
            args = docopt(__doc__, help=False, version=__version__)
        except DocoptExit:
            args = {"<command>": sys.argv[1]}
    else:
        print(__doc__)
    if args["<command>"]:
        args.update({"<tool>": os.path.basename(sys.argv[0])})
        # load <command> from entry_points
        for entry_point in working_set.iter_entry_points("%s.subcmd" % args["<tool>"]):
            if entry_point.name == args["<command>"]:
                subcommand = entry_point.load()
                sys.exit(subcommand())
        LOGGER.error(
            "'%s %s' is not a valid command", args["<tool>"], args["<command>"]
        )
        sys.exit(1)
    raise DocoptExit
