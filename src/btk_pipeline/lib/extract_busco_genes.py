#!/usr/bin/env python3
"""
Extract BUSCO gene sequences and set headers.

Usage: btk pipeline extract-busco-genes --busco PATH... --out FASTA

Options:
    --busco PATH         BUSCO full summary tsv file.
    --out FASTA          output FASTA filename.
"""

import codecs
import logging
import re
import sys
import tarfile
from pathlib import Path

from docopt import DocoptExit
from docopt import docopt
from tolkein import tofile

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
#         args["--busco"] = snakemake.input.busco
#         args["--out"] = snakemake.output.fasta
#     except NameError as err:
#         pass
#     for key, value in args.items:
#         sys.argv.append(key)
#         sys.argv.append(value)


def main():
    """Entry point."""
    try:
        # parse_args()
        args = docopt(__doc__)
    except DocoptExit:
        raise DocoptExit
    file_pattern = re.compile(r"busco_sequences\/(\w+?)_\w+\/(\d+at\d+).faa")
    header_pattern = re.compile(r">\S+:\d+-\d+")
    utf8reader = codecs.getreader("utf-8")
    try:
        with open(args["--out"], "w") as ofh:
            for busco_file in args["--busco"]:
                busco_dir = Path(busco_file).absolute().parent
                busco_seqs = "%s/busco_sequences.tar.gz" % busco_dir
                if not Path(busco_seqs).is_file():
                    continue
                tar = tarfile.open(busco_seqs)
                for tarinfo in tar.getmembers():
                    if tarinfo.name.endswith(".faa"):
                        match = file_pattern.match(tarinfo.name)
                        status, busco_id = match.groups()
                        with utf8reader(
                            tar.extractfile(
                                tarinfo,
                            )
                        ) as fh:
                            for line in fh:
                                if line.startswith(">"):
                                    header = line.strip()
                                    title = ""
                                    if header_pattern.match(header):
                                        title = "%s=%s=%s" % (header, busco_id, status)
                                    else:
                                        parts = header.split(" # ")
                                        title = "%s:%s-%s=%s=%s" % (
                                            re.sub(r"_[^_]+$", "", parts[0]),
                                            parts[1],
                                            parts[2],
                                            busco_id,
                                            status,
                                        )
                                    ofh.write("%s\n" % title)
                                else:
                                    ofh.write(line)
    except Exception as err:
        logger.error(err)
        exit(1)


if __name__ == "__main__":
    main()
