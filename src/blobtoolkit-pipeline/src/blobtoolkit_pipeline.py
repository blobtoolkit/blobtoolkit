#!/usr/bin/env python3

"""
BlobToolKit Pipeline.

Usage: blobtoolkit-pipeline [<command>] [<args>...] [-h|--help] [-v|--version]

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
See 'blobtoolkit-pipeline <command> --help' for more information on a specific command.

"""

import sys
from importlib.metadata import entry_points

from docopt import DocoptExit
from docopt import docopt
from lib.version import __version__


def iter_entry_points(group):
    """Return entry points for a group across supported Python versions."""
    discovered = entry_points()
    if hasattr(discovered, "select"):
        return discovered.select(group=group)
    return discovered.get(group, ())


def cli(rename=None):
    """Entry point."""
    # if sys.argv[1].match("-pipeline"):
    command = sys.argv[1]
    if len(sys.argv) > 2 and not sys.argv[2].startswith("-"):
        command = sys.argv[2]
    docs = __doc__
    if rename is not None:
        docs = docs.replace("blobtoolkit-pipeline", rename)
    try:
        args = docopt(docs, help=False, version=__version__)
    except DocoptExit:
        args = {"<command>": command}
    if args["<command>"]:
        # load <command> from entry_points
        for entry_point in iter_entry_points("blobtoolkit_pipeline.subcmd"):
            if entry_point.name == args["<command>"]:
                subcommand = entry_point.load()
                sys.exit(subcommand(rename))
    elif "--help" in args and args["--help"]:
        print(__doc__)
        exit(0)
    print(docs)
    raise DocoptExit
