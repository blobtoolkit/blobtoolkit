#!/usr/bin/env python3

"""
Transfer completed BlobDirs and intermediate files.

Usage:
  btk pipeline resume-pipeline --in PATH --out PATH

Options:
  --in PATH              Path to input directory.
  --out PATH             Path to output directory.
"""

import glob
import gzip
import re
import shutil
import tarfile
import zlib
from pathlib import Path
from time import sleep

import yaml
from docopt import docopt
from tolkein import tofile
from tolkein import tolog

LOGGER = tolog.logger(__name__)


def stream_gzip_decompress(stream):
    dec = zlib.decompressobj(32 + zlib.MAX_WBITS)  # offset 32 to skip the header
    for chunk in stream:
        rv = dec.decompress(chunk)
        if rv:
            yield rv


def untar_directory(archive, destdir, *, compress=False):
    LOGGER.info("Extracting files from %s to %s", archive, destdir)
    # make sure destdir exists
    Path(destdir).mkdir(parents=True, exist_ok=True)
    # extract files
    with tarfile.open(archive) as tfh:
        for member in tfh.getmembers():
            if ".stats." in member.name:
                f = tfh.extractfile(member)
                name = re.sub(r"^.+\/", "", member.name)
                name = name.replace(".gz", "")
                with open("%s/%s" % (destdir, name), "wb") as ofh:
                    for chunk in stream_gzip_decompress(f):
                        ofh.write(chunk)
    tools = [
        "windowmasker",
        "chunk_stats",
        "busco",
        "minimap",
        "cov_stats",
        "diamond",
        "diamond_blastp",
        "blastn",
        "window_stats",
        "blobtools",
        "view",
    ]
    for tool in tools:
        stats_file = Path("%s/%s.stats" % (destdir, tool))
        if stats_file.is_file():
            stats_file.touch()
            sleep(1)


def main():
    """Entry point."""
    opts = docopt(__doc__)
    indir = opts["--in"]
    config_file = "%s/config.yaml" % indir
    data = tofile.read_file(config_file)
    config = yaml.load(data, Loader=yaml.BaseLoader)
    destdir = opts["--out"]
    accession = config["assembly"]["accession"]
    Path("%s/%s" % (destdir, accession)).mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_file, "%s/%s/config.yaml" % (destdir, accession))
    tools = ["busco", "diamond", "diamond_blastp", "blastn", "window_stats"]
    for tool in tools:
        Path("%s/%s/%s" % (destdir, accession, tool)).mkdir(parents=True, exist_ok=True)
        for gzfile in glob.glob("%s/*%s*.gz" % (indir, tool)):
            with gzip.open(gzfile, "rb") as fh:
                file = gzfile.replace(".gz", "")
                file = file.replace(
                    "diamond_blastp.busco_genes.out", "diamond.busco_genes.out"
                )
                file = file.replace(
                    "diamond_blastx.reference_proteomes.out",
                    "diamond.reference_proteomes.out",
                )
                file = file.replace(indir, "%s/%s/%s" % (destdir, accession, tool))
                with open(file, "wb") as ofh:
                    shutil.copyfileobj(fh, ofh)

    for tf in glob.glob("%s/*.pipeline.tar" % indir):
        untar_directory(tf, "%s/%s" % (destdir, accession))
    for tf in glob.glob("%s/*.busco*.tar" % indir):
        with tarfile.open(tf) as tfh:
            tfh.extractall(path="%s/%s/busco" % (destdir, accession))
    # print(data)


if __name__ == "__main__":
    main()
