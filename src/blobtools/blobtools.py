#!/usr/bin/env python3

"""
BlobTools2 - assembly exploration, QC and filtering.

usage: blobtools [<command>] [<args>...] [-h|--help] [--version]

commands:
    add             add data to a BlobDir
    create          create a new BlobDir
    filter          filter a BlobDir
    host            host interactive view of all BlobDirs in a directory
    replace         call blobtools add with --replace flag
    remove          remove one or more fields from a BlobDir
    validate        validate a BlobDir
    view            generate plots using BlobToolKit Viewer
    -h, --help      show this
    -v, --version   show version number
See 'blobtools <command> --help' for more information on a specific command.

examples:
    # 1. Create a new BlobDir from a FASTA file:
    blobtools create --fasta examples/assembly.fasta BlobDir

    # 2. Create a new BlobDir from a BlobDB:
    blobtools create --blobdb examples/blobDB.json BlobDir

    # 3. Add Coverage data from a BAM file:
    blobtools add --cov examples/assembly.reads.bam BlobDir

    # 4. Assign taxonomy from BLAST hits:
    blobtools add add --hits examples/blast.out --taxdump ../taxdump BlobDir

    # 5. Add BUSCO results:
    blobtools add --busco examples/busco.tsv BlobDir

    # 6. Host an interactive viewer:
    blobtools host BlobDir

    # 7. Filter a BlobDir:
    blobtools filter --param length--Min=5000 --output BlobDir_len_gt_5000 BlobDir
"""


import os
import sys
import warnings
from datetime import datetime
from importlib.metadata import entry_points

from docopt import DocoptExit
from docopt import docopt

from .lib.tolkein_compat import tolog
from .lib.version import __version__

LOGGER = tolog.logger(__name__)

# Python EOL dates from https://devguide.python.org/versions/
PYTHON_EOL_DATES = {
    (3, 9): datetime(2025, 10, 5),
    (3, 10): datetime(2026, 10, 4),
    (3, 11): datetime(2027, 10, 24),
    (3, 12): datetime(2028, 10, 2),
    (3, 13): datetime(2029, 10, 1),
}

# Minimum supported Python version
MIN_PYTHON_VERSION = (3, 10)


def check_python_version():
    """Check for deprecated or EOL Python versions and warn users."""
    current_version = (sys.version_info.major, sys.version_info.minor)
    current_date = datetime.now()

    # Check if running on unsupported version
    if current_version < MIN_PYTHON_VERSION:
        warnings.warn(
            f"Python {current_version[0]}.{current_version[1]} is no longer supported by BlobToolKit. "
            f"Please upgrade to Python {MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]} or later.",
            DeprecationWarning,
            stacklevel=3,
        )
        return

    # Check if running on version approaching EOL
    if current_version in PYTHON_EOL_DATES:
        eol_date = PYTHON_EOL_DATES[current_version]
        days_until_eol = (eol_date - current_date).days

        if days_until_eol < 0:
            # Already past EOL
            warnings.warn(
                f"Python {current_version[0]}.{current_version[1]} reached end-of-life on {eol_date.strftime('%Y-%m-%d')}. "
                "BlobToolKit support will be removed in a future release. "
                f"Please upgrade to Python {MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]} or later for continued support and security updates.",
                DeprecationWarning,
                stacklevel=3,
            )
        elif days_until_eol < 180:  # Warn 6 months before EOL
            warnings.warn(
                f"Python {current_version[0]}.{current_version[1]} will reach end-of-life on {eol_date.strftime('%Y-%m-%d')} "
                f"({days_until_eol} days). "
                f"Please plan to upgrade to Python {MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]} or later.",
                FutureWarning,
                stacklevel=3,
            )


check_python_version()


def iter_entry_points(group):
    """Yield entry points for a given group across Python versions."""
    eps = entry_points()
    if hasattr(eps, "select"):
        return eps.select(group=group)
    return eps.get(group, [])


def suggest_option(command):
    options = {
        "host": ["host"],
    }
    options_list = options.get(command, [])
    options_list.append("full")
    LOGGER.error(
        f"The `blobtools {command}` command is not available, to enable this option use{' one of' if len(options_list) > 1 else ''}:"
    )
    for option in options_list:
        print(
            f"                                - pip install blobtoolkit[{option}]",
            file=sys.stderr,
        )


def cli():
    """Entry point."""
    args = {}
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
    if "<command>" in args and args["<command>"]:
        args.update({"<tool>": os.path.basename(sys.argv[0])})
        # load <command> from entry_points
        command_index = sys.argv.index(args["<command>"])
        if args["<command>"] in ["create", "replace"]:
            if len(sys.argv) > command_index + 1:
                if args["<command>"] == "replace" and "--replace" not in sys.argv:
                    sys.argv.insert(command_index + 1, "--replace")
                args["<command>"] = "add"
            sys.argv[command_index] = "add"
        for entry_point in iter_entry_points(f"{args['<tool>']}.subcmd"):
            if entry_point.name == args["<command>"]:
                try:
                    subcommand = entry_point.load()
                    sys.exit(subcommand())
                except ModuleNotFoundError:
                    suggest_option(args["<command>"])
                    exit(1)
        LOGGER.error(
            "'%s %s' is not a valid command", args["<tool>"], args["<command>"]
        )
        sys.exit(1)
    elif "--help" in args and args["--help"]:
        print(__doc__)
        exit(0)
    elif "--version" in args and args["--version"]:
        print(__version__)
        exit(0)
    raise DocoptExit
