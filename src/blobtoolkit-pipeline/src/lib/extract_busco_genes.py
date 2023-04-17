#!/usr/bin/env python3
"""
Extract BUSCO gene sequences and set headers.

Usage: blobtoolkit-pipeline extract-busco-genes --busco PATH... --out FASTA

Options:
    --busco PATH         BUSCO full summary tsv file or busco_sequences directory.
    --out FASTA          output FASTA filename.
"""

import codecs
import glob
import logging
import re
import tarfile
from pathlib import Path

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
#         args["--busco"] = snakemake.input.busco
#         args["--out"] = snakemake.output.fasta
#     except NameError as err:
#         pass
#     for key, value in args.items:
#         sys.argv.append(key)
#         sys.argv.append(value)


def parse_busco_file(fh, ofh, status, busco_id):
    """Parse a busco sequences file."""
    header_pattern = re.compile(r">\S+:\d+-\d+")
    for line in fh:
        if line.startswith(">"):
            header = line.strip()
            title = ""
            if header_pattern.match(header):
                title = f"{header}={busco_id}={status}"
            else:
                parts = header.split(" # ")
                title = f'{re.sub("_[^_]+$", "", parts[0])}:{parts[1]}-{parts[2]}={busco_id}={status}'
            ofh.write("%s\n" % title)
        else:
            ofh.write(line)


def main(rename=None):
    """Entry point."""
    docs = __doc__
    if rename is not None:
        docs = docs.replace("blobtoolkit-pipeline", rename)
    try:
        args = docopt(docs)
    except DocoptExit as e:
        raise DocoptExit from e
    file_pattern = re.compile(r"busco_sequences\/(\w+?)_\w+\/(\d+at\d+)\.faa")
    try:
        with open(args["--out"], "w") as ofh:
            for busco_file in args["--busco"]:
                if busco_file.endswith("busco_sequences"):
                    # plain directory
                    parse_busco_directory(file_pattern, ofh, busco_file)
                else:
                    # tarred archive
                    utf8reader = codecs.getreader("utf-8")
                    busco_dir = Path(busco_file).absolute().parent
                    busco_seqs = f"{busco_dir}/busco_sequences.tar.gz"
                    if not Path(busco_seqs).is_file():
                        continue
                    parse_tarred_busco_directory(
                        file_pattern, ofh, utf8reader, busco_seqs
                    )

    except Exception as err:
        logger.error(err)
        exit(1)


def parse_tarred_busco_directory(file_pattern, ofh, utf8reader, busco_seqs):
    """Parse a tarred busco directory."""
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
                parse_busco_file(fh, ofh, status, busco_id)


def parse_busco_directory(file_pattern, ofh, busco_file):
    """Parse a plain busco_directory."""
    files = glob.glob(f"{busco_file}/*/*.faa")
    for file in files:
        match = file_pattern.search(file)
        status, busco_id = match.groups()
        with open(file, "r") as fh:
            parse_busco_file(fh, ofh, status, busco_id)


if __name__ == "__main__":
    main()
