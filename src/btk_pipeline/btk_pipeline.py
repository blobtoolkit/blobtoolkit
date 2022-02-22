#!/usr/bin/env python3

"""
BlobToolKit Pipeline.

Usage: blobtools pipeline [<command>] [<args>...] [-h|--help] [--version]

Commands:
    data                      Fetch pipeline data
    run                       Run BlobToolKit Pipeline
    add-summary-to-metadata   Pipeline helper script
    chunk-fasta               Pipeline helper script
    count-busco-genes         Pipeline helper script
    extract-busco-genes       Pipeline helper script
    generate-config           Pipeline helper script
    generate-static-images    Pipeline helper script
    transfer-completed        Pipeline helper script
    unchunk-blast             Pipeline helper script
    window-stats              Pipeline helper script
    -h, --help      Show this
    -v, --version   Show version number
See 'blobtools pipeline <command> --help' for more information on a specific command.

"""

import sys

from docopt import DocoptExit
from docopt import docopt
from pkg_resources import working_set

from .lib.version import __version__


def main():
    """Entry point."""
    if len(sys.argv) > 2:
        try:
            args = docopt(__doc__, help=False, version=__version__)
        except DocoptExit:
            args = {"<command>": sys.argv[2]}
        if args["<command>"]:
            # load <command> from entry_points
            for entry_point in working_set.iter_entry_points("pipeline.subcmd"):
                if entry_point.name == args["<command>"]:
                    subcommand = entry_point.load()
                    sys.exit(subcommand())
    print(__doc__)
    raise DocoptExit
