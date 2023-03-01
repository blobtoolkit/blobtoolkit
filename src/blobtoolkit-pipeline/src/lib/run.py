#!/usr/bin/env python3
"""
Run BlobToolKit Pipeline.

Usage: blobtoolkit-pipeline run [--dry-run] [--tool STRING] [--unlock]
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
    snakefile = os.path.join(script_dir, "blobtoolkit-snakefiles", "blobtoolkit.smk")
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
    return exit_code


def run_pipeline(workdir, args):
    """Run Snakemake pipeline."""
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    snakefile = os.path.join(
        script_dir, "blobtoolkit-snakefiles", f'{args["--tool"]}.smk'
    )
    dry_run = "-n" if "--dry-run" in args and args["--dry-run"] is True else ""
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
    return exit_code


def run_snakemake_pipeline(args):
    config = os.path.abspath(args["--config"])
    workdir = (
        os.path.abspath(args["--workdir"])
        if args["--workdir"] is not None
        else os.path.abspath(os.path.dirname(config))
    )
    if config != f"{workdir}/config.yaml":
        shutil.copy2(config, f"{workdir}/config.yaml")
    if args["--unlock"]:
        exit_code = unlock_working_directory(workdir)
        exit(exit_code)
    exit_code = run_pipeline(workdir, args)
    exit(exit_code)


def main(rename=None):
    """Entry point."""
    docs = __doc__
    if rename is not None:
        docs = docs.replace("blobtoolkit-pipeline", rename)
    try:
        args = docopt(docs)
    except DocoptExit as e:
        raise DocoptExit from e
    try:
        run_snakemake_pipeline(args)
    except Exception as err:
        logger.error(err)
        exit(1)
