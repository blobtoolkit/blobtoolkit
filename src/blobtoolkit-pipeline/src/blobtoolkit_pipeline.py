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
import warnings
from datetime import datetime
from importlib.metadata import entry_points

from docopt import DocoptExit
from docopt import docopt
from lib.version import __version__

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
