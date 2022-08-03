#!/usr/bin/env python3
"""
Run BlobToolKit Pipeline.

Usage: btk pipeline run [--dry-run] [--tool STRING] [--unlock]
                        [--workdir DIR] [--threads INT] --config YAML

Options:
    --config YAML  YAML format configuration filename.
    --dry-run      Dry run flag.
    --threads INT  Number of threads to use. [Default: 32]
    --tool STRING  Pipeline tool to run. [Default: blobtoolkit]
    --unlock       Flag to unlock working directory.
    --workdir DIR  Full Path of working directory.
                   Default is the directory containing config file.
"""

import logging
import os
import shlex
import shutil
import signal
import subprocess

from docopt import DocoptExit
from docopt import docopt

logger_config = {
    "level": logging.INFO,
    "format": "%(asctime)s [%(levelname)s] line %(lineno)d %(message)s",
    "filemode": "w",
}
logging.basicConfig(**logger_config)
logger = logging.getLogger()


def run_command(cmd):
    """Run a command and capture CTRL-C."""
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    p = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, cwd=script_dir)
    exit_code = p.wait()
    return exit_code


def unlock_working_directory(workdir):
    """Unlock Snakemake working directory after a failed run."""
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    snakefile = os.path.join(script_dir, "blobtoolkit.smk")
    cmd = """
    snakemake -p \
          -j 1 \
          --directory %s/blobtoolkit \
          --configfile %s/config.yaml \
          --unlock \
          -s %s
    """ % (
        workdir,
        workdir,
        snakefile,
    )
    exit_code = run_command(cmd)


def run_pipeline(workdir, args):
    """Run Snakemake pipeline."""
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    snakefile = os.path.join(script_dir, "%s.smk" % args["--tool"])
    if "--dry-run" in args and args["--dry-run"] is True:
        dry_run = "-n"
    else:
        dry_run = ""
    cmd = """
    snakemake -p %s \
          -j %s \
          --directory %s/blobtoolkit \
          --configfile %s/config.yaml \
          -s %s
    """ % (
        dry_run,
        args["--threads"],
        workdir,
        workdir,
        snakefile,
    )
    exit_code = run_command(cmd)


def main():
    """Entry point."""
    try:
        args = docopt(__doc__)
    except DocoptExit:
        raise DocoptExit
    try:
        config = os.path.abspath(args["--config"])
        if args["--workdir"] is not None:
            workdir = os.path.abspath(args["--workdir"])
        else:
            workdir = os.path.abspath(os.path.dirname(config))
        if config != "%s/config.yaml" % workdir:
            shutil.copy2(config, "%s/config.yaml" % workdir)
        if args["--unlock"]:
            exit_code = unlock_working_directory(workdir)
            exit(exit_code)
        exit_code = run_pipeline(workdir, args)
        exit(exit_code)

    except Exception as err:
        logger.error(err)
        exit(1)
