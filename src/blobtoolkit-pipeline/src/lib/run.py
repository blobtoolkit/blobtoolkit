#!/usr/bin/env python3
"""
Run BlobToolKit Pipeline.

Usage: blobtoolkit-pipeline run [options]

Options:
    --config YAML  YAML config file.  In Snakemake mode: path to the pipeline
                   config.yaml.  In --nextflow mode: optional YAML file that
                   can supply any of the Nextflow options listed below; CLI
                   flags always override file values.  Supported flat keys:
                   accession, fasta, taxon, input, outdir, taxdump, blastp,
                   blastn, blastx, busco_lineages, profile, revision, workflow,
                   align, nextflow_binary, dry_run.  Nested Snakemake keys
                   assembly.accession, taxon.taxid, and busco.lineages are
                   also recognised.
    --dry-run      Dry run flag.
    --threads INT  Number of threads to use. [Default: 32]
    --tool STRING  Pipeline tool to run. [Default: blobtoolkit]
    --unlock       Flag to unlock working directory.
    --workdir DIR  Full Path of working directory.
                   Default is the directory containing config file.
    --nextflow     Run the preferred external Nextflow workflow.
    --input CSV    Sample sheet in CSV format for Nextflow workflow.
    --fasta FASTA  Assembly FASTA path for Nextflow workflow.
    --accession STRING  Assembly accession for Nextflow workflow.
    --taxon STRING      Taxon identifier for Nextflow workflow.
    --outdir DIR        Output directory for Nextflow workflow.
    --taxdump PATH      Path to NCBI taxdump directory.
    --blastp PATH       Path to blastp diamond database.
    --blastn PATH       Path to blastn database index.
    --blastx PATH       Path to blastx diamond database.
    --profile STRING    Nextflow profile(s). [Default: sanger,singularity]
    --revision STRING   Nextflow workflow revision. [Default: 0.10.0]
    --workflow STRING   Nextflow workflow name. [Default: sanger-tol/blobtoolkit]
    --nextflow-binary STRING  Nextflow executable to invoke. [Default: nextflow]
    --align        Enable alignment stage in Nextflow workflow.
    --busco-lineages STRING  Comma-separated BUSCO lineage(s) to use, e.g.
                             nematoda_odb10,eukaryota_odb10.
"""

import logging
import os
import shlex
import shutil
import subprocess
import warnings

import yaml
from docopt import DocoptExit
from docopt import docopt

# ---------------------------------------------------------------------------
# Mapping from YAML config keys -> docopt arg names.
# Supports both a flat BTK-Nextflow style and the nested Snakemake style so
# that existing config.yaml files can be reused with minimal changes.
# ---------------------------------------------------------------------------
_CONFIG_KEY_MAP = {
    # flat keys (new / preferred)
    "accession": "--accession",
    "fasta": "--fasta",
    "taxon": "--taxon",
    "input": "--input",
    "outdir": "--outdir",
    "taxdump": "--taxdump",
    "blastp": "--blastp",
    "blastn": "--blastn",
    "blastx": "--blastx",
    "profile": "--profile",
    "revision": "--revision",
    "workflow": "--workflow",
    "align": "--align",
    "busco_lineages": "--busco-lineages",
    "nextflow_binary": "--nextflow-binary",
    "dry_run": "--dry-run",
}

# Nested Snakemake config keys expressed as dot-paths → docopt arg names.
_NESTED_KEY_MAP = {
    "assembly.accession": "--accession",
    "taxon.taxid": "--taxon",
    "assembly.fasta": "--fasta",
    "busco.lineages": "--busco-lineages",
}


def _get_nested(data, dotpath):
    """Return value at a dot-separated path in a nested dict, or None."""
    parts = dotpath.split(".")
    node = data
    for part in parts:
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


def load_nextflow_config(config_path):
    """Load a YAML config file and return a dict of {--flag: value}.

    Supports both a flat BTK-Nextflow style config and the nested Snakemake
    config.yaml format.  Flat keys take precedence over nested mappings.
    """
    with open(config_path) as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"Config file {config_path!r} must contain a YAML mapping.")

    result = {}

    # Nested Snakemake-style mappings first (lower priority).
    for dotpath, flag in _NESTED_KEY_MAP.items():
        val = _get_nested(data, dotpath)
        if val is not None:
            result[flag] = ",".join(val) if isinstance(val, list) else str(val)

    # Flat keys (higher priority, overwrite nested values).
    for key, flag in _CONFIG_KEY_MAP.items():
        if key in data and data[key] is not None:
            val = data[key]
            # Boolean flags (--align, --dry-run) become True/False.
            result[flag] = True if isinstance(val, bool) else str(val)

    return result


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
    return p.wait()


def run_list_command(cmd, cwd=None):
    """Run a tokenized command and capture CTRL-C."""
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, cwd=cwd)
    return p.wait()


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
    return run_command(cmd)


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
    return run_command(cmd)


def run_snakemake_pipeline(args):
    if not args["--config"]:
        raise ValueError("--config YAML is required for Snakemake mode.")
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


def run_nextflow_pipeline(args):
    """Run external Nextflow BlobToolKit workflow."""
    # Merge config file values first; CLI args take precedence.
    if args.get("--config"):
        try:
            file_args = load_nextflow_config(args["--config"])
        except Exception as exc:
            raise ValueError(f"Failed to load config file: {exc}") from exc
        for flag, value in file_args.items():
            # Only apply if the CLI did not supply an explicit value.
            if not args.get(flag):
                args[flag] = value

    required = ["--input", "--fasta", "--accession", "--taxon", "--outdir"]
    if missing := [flag for flag in required if not args.get(flag)]:
        raise ValueError(
            "The following required arguments are missing for --nextflow mode: "
            + ", ".join(missing)
        )

    nextflow_binary = args["--nextflow-binary"]
    if shutil.which(nextflow_binary) is None:
        raise FileNotFoundError(
            f"Unable to find Nextflow executable '{nextflow_binary}' in PATH"
        )

    cmd = [
        nextflow_binary,
        "run",
        args["--workflow"],
        "-r",
        args["--revision"],
        "-profile",
        args["--profile"],
        "--input",
        os.path.abspath(args["--input"]),
        "--fasta",
        os.path.abspath(args["--fasta"]),
        "--accession",
        args["--accession"],
        "--taxon",
        args["--taxon"],
        "--outdir",
        os.path.abspath(args["--outdir"]),
    ]

    optional_path_flags = ["--taxdump", "--blastp", "--blastn", "--blastx"]
    for flag in optional_path_flags:
        if value := args.get(flag):
            cmd.extend([flag, os.path.abspath(value)])

    if value := args.get("--busco-lineages"):
        cmd.extend(["--busco_lineages", value])

    if args["--align"]:
        cmd.append("--align")
    if args["--dry-run"]:
        cmd.append("-stub-run")
    if args["--unlock"]:
        warnings.warn(
            "--unlock is not supported for Nextflow runs and will be ignored.",
            UserWarning,
            stacklevel=2,
        )
    if args["--workdir"]:
        warnings.warn(
            "--workdir is not used by Nextflow handoff and will be ignored.",
            UserWarning,
            stacklevel=2,
        )
    if args["--threads"] != "32":
        warnings.warn(
            "--threads is not currently mapped to Nextflow resources and will be ignored.",
            UserWarning,
            stacklevel=2,
        )

    logger.info("Running Nextflow command: %s", " ".join(shlex.quote(t) for t in cmd))
    exit_code = run_list_command(cmd)
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
        if args["--nextflow"]:
            run_nextflow_pipeline(args)
        # Warn about deprecated Snakemake pipeline only when using legacy mode.
        warnings.warn(
            "The Snakemake-based BlobToolKit pipeline is no longer actively maintained. "
            "Please use the actively supported Nextflow implementation instead: "
            "https://pipelines.tol.sanger.ac.uk/blobtoolkit. "
            "To run the Nextflow pipeline via this command, use the --nextflow flag and provide the required arguments.",
            FutureWarning,
            stacklevel=2,
        )
        run_snakemake_pipeline(args)
    except Exception as err:
        logger.error(err)
        exit(1)
