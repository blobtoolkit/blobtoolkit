#!/usr/bin/env python3
"""
Unchunk BLAST results.

Update coordinates and remove seq name suffix from chunked blast results.

Usage: blobtoolkit-pipeline unchunk-blast [--count INT] --in TSV --out TSV

Options:
    --in TSV     input filename.
    --out TSV    output filename.
    --count INT  number of results to keep per chunk. [Default: 10]
"""

import logging
import re
import sys
from collections import defaultdict

from docopt import DocoptExit
from docopt import docopt

logger_config = {
    "level": logging.INFO,
    "format": "%(asctime)s [%(levelname)s] line %(lineno)d %(message)s",
    "filemode": "w",
}
logging.basicConfig(**logger_config)
logger = logging.getLogger()


# def parse_args():
#     """Parse snakemake args if available."""
#     args = {}
#     try:
#         args["--in"] = snakemake.input[0]
#         args["--count"] = str(snakemake.params.max_target_seqs)
#         args["--out"] = snakemake.output[0]
#         for key, value in args.items:
#             sys.argv.append(key)
#             sys.argv.append(value)
#     except NameError as err:
#         pass


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
        lines = defaultdict(dict)
        chunk_counts = defaultdict(int)
        with open(args["--in"], "r") as fh:
            for line in fh.readlines():
                if "\n" not in line:
                    line += "\n"
                fields = line.split("\t")
                if fields[0]:
                    name, start = re.split("_-_", fields[0])
                    fields[0] = name
                    fields[3] = name
                    fields[9] = str(int(fields[9]) + int(start))
                    fields[10] = str(int(fields[10]) + int(start))
                    if start not in lines[name]:
                        lines[name][start] = []
                        chunk_counts[name] += 1
                    lines[name][start].append("\t".join(fields))
        with open(args["--out"], "w") as ofh:
            for name in chunk_counts.keys():
                length = len(lines[name])
                n = int(args["--count"])
                for start in lines[name].keys():
                    for i in range(n):
                        if i < len(lines[name][start]):
                            ofh.write("%s" % (lines[name][start][i]))
    except Exception as err:
        logger.error(err)
        exit(1)
