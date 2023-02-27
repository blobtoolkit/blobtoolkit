#!/usr/bin/env python3

"""
Add summary data to metadata.

Usage: blobtoolkit-pipeline add-summary-to-metadata --config FILE --out FILE

Options:
    --config FILE   YAML format config file
    --out FILE      Output file
"""


import logging
import re
import subprocess
import traceback
from collections import OrderedDict

import yaml
from docopt import DocoptExit
from docopt import docopt
from tolkein import tofile

from .functions import read_similarity_settings
from .functions import reads_by_prefix
from .version import __version__

logger_config = {
    "level": logging.INFO,
    "format": "%(asctime)s [%(levelname)s] line %(lineno)d %(message)s",
    "filemode": "w",
}
logging.basicConfig(**logger_config)
logger = logging.getLogger()

logger.info(f"Starting script: {__file__}")


def add_pipeline_version(meta):
    """Add pipeline version info to metadata."""
    meta["settings"]["pipeline"] = "https://github.com/blobtoolkit/blobtoolkit"
    meta["settings"]["release"] = __version__


def add_software_versions(meta):
    """Add software versions to meta."""
    programs = {
        "blastn": {"flag": "-version", "regex": r":\s*(\S+)$"},
        "blobtools": {"regex": r"\s*v(\S+)$"},
        "busco": {"regex": r"\s+(\S+)$"},
        "diamond": {"regex": r"version\s*(\S+)$"},
        "minimap2": {"regex": r"(\S+)"},
        "blobtk": {"regex": r"\s+(\S+)$"},
        "python": {"regex": r"\s+(\S+)$"},
        "samtools": {"regex": r"\s+(\S+)$"},
        "seqtk": {"flag": "", "line": 2, "status": 1, "regex": r":\s*(\S+)$"},
        "snakemake": {"regex": r"(\S+)"},
    }
    versions = {}
    logger.info("Adding software version details")
    for key, opts in programs.items():
        cmd = [key]
        if flag := opts.get("flag", "--version"):
            cmd.append(flag)
        try:
            p = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                encoding="utf-8",
            )
            (output, err) = p.communicate()
            status = p.wait()
            expected_status = opts.get("status", 0)
            if expected_status == 1:
                output = err
                status = 0
            if status == 0:
                for idx, line in enumerate(output.split("\n")):
                    if idx == opts.get("line", 0):
                        match = re.search(opts["regex"], line)
                        version = match.groups()[0]
            else:
                version = "unknown"
        except FileNotFoundError:
            version = "not found"
        logger.info("%s: %s", key, version)
        versions[key] = version
    meta["settings"]["software_versions"] = versions


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
        meta = {}
        config = args["--config"]
        if not isinstance(config, dict):
            config = tofile.load_yaml(config)
        meta["assembly"] = config["assembly"]
        meta["taxon"] = config["taxon"]
        meta["settings"] = {
            "blast_chunk": 100000,
            "blast_max_chunks": 10,
            **config["settings"],
        }
        meta["similarity"] = {
            "diamond_blastx": read_similarity_settings(config, "diamond_blastx"),
            "diamond_blastp": read_similarity_settings(config, "diamond_blastp"),
        }
        meta["reads"] = reads_by_prefix(config)
        for key, value in meta["reads"].items():
            if "url" in value and not isinstance(value["url"], list):
                value["url"] = value["url"].split(";")
            elif "url" not in value:
                value["url"] = []
        add_pipeline_version(meta)
        add_software_versions(meta)
        meta["version"] = config.get("version", 1)
        meta["revision"] = config.get("revision", 0)

        yaml.add_representer(
            OrderedDict,
            lambda dumper, data: dumper.represent_mapping(
                "tag:yaml.org,2002:map", data.items()
            ),
        )
        yaml.add_representer(
            tuple,
            lambda dumper, data: dumper.represent_sequence(
                "tag:yaml.org,2002:seq", data
            ),
        )
        yaml.Dumper.ignore_aliases = lambda *args: True
        with open(args["--out"], "w") as fh:
            fh.write(yaml.dump(meta, indent=1))

    except Exception as err:
        logger.error(traceback.format_exc())
        exit(13)


if __name__ == "__main__":
    main()
