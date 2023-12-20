#!/usr/bin/env python3
"""
Split long sequences into chunks.

Usage: blobtoolkit-pipeline chunk-fasta --in FASTA [--chunk INT] [--overlap INT] [--max-chunks INT]
                [--busco TSV] [--min-length INT] [--out FASTA] [--bed BEDFILE]

Options:
    --in FASTA           input FASTA file.
    --busco TSV          BUSCO full summary tsv file.
    --chunk INT          sequences greater than CHUNK bp will be split. [Default: 100000]
    --overlap INT        length of overlap when splitting sequences. [Default: 0]
    --max-chunks INT     maximum number of chunks to split a sequence into. [Default: 10]
    --min-length INT     minimum sequence length. [Default: 1000]
    --out FASTA          output FASTA filename or suffix. [Default: .chunked]
    --bed BEDFILE        output BED/BED-like TSV filename or suffix. [Default: .bed]
"""

import logging
import math
import re
import shlex
import sys
from collections import Counter
from collections import defaultdict
from itertools import groupby
from subprocess import PIPE
from subprocess import Popen

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


def chunk_size(value):
    """Calculate nice value for chunk size."""
    return round(value / 1000 + 0.5) * 1000


def chunk_fasta(
    fastafile, *, chunk=math.inf, overlap=0, max_chunks=math.inf, min_length=0
):
    """Read FASTA file one sequence at a time and split long sequences into chunks."""
    cmd = "cat %s" % fastafile
    if fastafile.endswith(".gz"):
        cmd = "zcat %s" % fastafile
    title = ""
    seq = ""
    segment = chunk + overlap
    with Popen(shlex.split(cmd), encoding="utf-8", stdout=PIPE, bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == ">"))
        for header in faiter:
            title = header.__next__()[1:].strip().split()[0]
            seq = "".join(map(lambda s: s.strip(), faiter.__next__()))
            seq_length = len(seq)
            my_chunk = chunk
            if seq_length > segment:
                n = (seq_length + chunk) // chunk
                if n > max_chunks:
                    my_chunk = chunk_size(seq_length / max_chunks)
                    n = max_chunks
                for i in range(0, seq_length, my_chunk):
                    subseq = seq[i : i + my_chunk + overlap]
                    yield {
                        "title": title,
                        "seq": subseq,
                        "chunks": n,
                        "start": i,
                        "end": i + my_chunk + overlap,
                        "length": my_chunk + overlap,
                    }
            elif seq_length > min_length:
                yield {
                    "title": title,
                    "seq": seq,
                    "chunks": 1,
                    "start": 0,
                    "end": seq_length,
                    "length": seq_length,
                }


def check_for_unmasked_bases(seq, min_unmasked=20):
    """Check sequence has runs of unmasked bases."""
    return bool(re.search(r"[ACGT]{" + str(min_unmasked) + "}", seq))


def check_for_masked_bases(seq, max_masked=20):
    """Check sequence has runs of masked bases."""
    return bool(re.search(r"[acgtnN]{" + str(max_masked) + "}", seq))


def parse_busco_full_summary(busco_file, chunk=100000):
    """Parse a busco full summary file."""
    logger.info("Parsing BUSCO full summary file")
    locations = defaultdict(list)
    with tofile.open_file_handle(busco_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[1] in ("Complete", "Duplicated"):
                if ":" in parts[2] and parts[5] == "-":
                    start = int(parts[4])
                    end = int(parts[3])
                else:
                    start = int(parts[3])
                    end = int(parts[4])
                locations[parts[2].split(":")[0]].append((start, end))
    windows = {}
    for title, tuples in locations.items():
        tuples.sort(key=lambda tup: tup[0])
        windows[title] = []
        start_index = 0
        for location in tuples:
            windows[title].append([location[0], location[0] + chunk, 0])
            for window in windows[title][start_index:]:
                if location[1] < window[1]:
                    window[2] += 1
                else:
                    start_index += 1
        windows[title].sort(key=lambda window: window[2], reverse=True)
    return windows


def chunk_by_busco(seq, seqs, busco_windows, args):
    """Use BUSCO gene density to define chunks."""
    subseq_found = False
    has_busco = False
    if seq["title"] in busco_windows:
        for window in busco_windows[seq["title"]]:
            if seq["start"] <= window[0] and seq["end"] >= window[1]:
                has_busco = True
                subseq = seq["seq"][window[0] - seq["start"] : window[1] - seq["start"]]
                # check there are runs of unmasked bases
                if check_for_unmasked_bases(subseq, 1000):
                    seq["seq"] = subseq
                    seq["start"] = window[0]
                    subseq_found = True
                    break
    if not has_busco:
        # extract subseq from middle of chunk
        chunk = int(args["--chunk"])
        offset = int((seq["length"] - chunk) / 2)
        if offset > 0:
            midseq = seq["seq"][offset:-offset]
            # check there are runs of unmasked bases
            if check_for_unmasked_bases(midseq, 1000):
                seq["seq"] = midseq
                seq["start"] += offset
                seq["end"] -= offset
                subseq_found = True
            else:
                # walk along sequence to find a chunk with fewer masked bases
                while offset > chunk:
                    offset -= chunk
                    lseq = seq["seq"][offset : offset + chunk]
                    if check_for_unmasked_bases(lseq, 1000):
                        seq["seq"] = lseq
                        seq["start"] += offset
                        seq["end"] = seq["start"] + chunk
                        subseq_found = True
                        break
                    rseq = seq["seq"][-offset - chunk : -offset]
                    if check_for_unmasked_bases(rseq, 1000):
                        seq["seq"] = rseq
                        seq["start"] = seq["end"] - offset - chunk
                        seq["end"] -= offset
                        subseq_found = True
                        break
    if subseq_found:
        seqs.append((seq))


def seq_stats(seq):
    """Calculate basic sequence statistics."""
    counter = Counter(seq)
    gc = 0
    at = 0
    n = 0
    masked = 0
    for base in ["g", "c", "G", "C"]:
        if base in counter:
            gc += counter[base]
    for base in ["a", "t", "A", "T"]:
        if base in counter:
            at += counter[base]
    for base in ["n", "N"]:
        if base in counter:
            n += counter[base]
    for base in ["a", "c", "g", "t", "n"]:
        if base in counter:
            masked += counter[base]
    atgc = at + gc
    if atgc > 0:
        gc = gc / atgc
    else:
        gc = 0
    return {"gc": gc, "n": n / len(seq), "ncount": n, "masked": masked / len(seq)}


def write_bedfiles(bed_data, args):
    """Write bedfiles."""
    if "--bed" not in args or not args["--bed"]:
        return
    locations = []
    rows = []
    header = ""
    for title, arr in bed_data.items():
        for obj in arr:
            if not header:
                header = "sequence\tstart\tend\t%s\n" % "\t".join(obj["stats"].keys())
                rows.append(header)
            row = "%s\t%d\t%d" % (title, obj["start"], obj["end"])
            locations.append("%s\n" % row)
            for key, value in obj["stats"].items():
                if key != "ncount":
                    row += "\t%.4f" % value
                else:
                    row += "\t%d" % value
            rows.append("%s\n" % row)
    bedfile = args["--bed"]
    if bedfile.startswith("."):
        bedfile = "%s%s" % (args["--in"], bedfile)
    bedfile = re.sub(r"\.bed$", "", bedfile)
    bedfile = re.sub(r"\.tsv$", "", bedfile)
    with open("%s.mask.bed" % bedfile, "w") as ofh:
        ofh.writelines(locations)
    with open("%s.tsv" % bedfile, "w") as ofh:
        ofh.writelines(rows)


# def parse_args():
#     """Parse snakemake args if available."""
#     args = {}
#     try:
#         args["--in"] = snakemake.input.fasta
#         args["--chunk"] = str(snakemake.params.chunk)
#         args["--overlap"] = str(snakemake.params.overlap)
#         args["--max-chunks"] = str(snakemake.params.max_chunks)
#         args["--min-length"] = str(snakemake.params.min_length)
#         try:
#             args["--busco"] = snakemake.input.busco
#         except AttributeError:
#             args["--busco"] = "None"
#         try:
#             args["--out"] = snakemake.output.fasta
#         except AttributeError:
#             args["--out"] = "None"
#         try:
#             args["--bed"] = snakemake.params.bed
#         except AttributeError:
#             args["--bed"] = "None"
#         for key, value in args.items():
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
        make_chunks(args)
    except Exception as err:
        logger.error(err)
        exit(1)


def make_chunks(args):
    busco_windows = (
        parse_busco_full_summary(args["--busco"])
        if (
            "--busco" in args
            and args["--busco"] is not None
            and args["--busco"] != "None"
        )
        else {}
    )
    logger.info(f'Splitting {args["--in"]} into chunks')
    bed_data = defaultdict(list)
    seqs = []
    for seq in chunk_fasta(
        args["--in"],
        chunk=int(args["--chunk"]),
        overlap=int(args["--overlap"]),
        max_chunks=int(args["--max-chunks"]),
    ):
        if (
            busco_windows
            and seq["chunks"] == int(args["--max-chunks"])
            and seq["length"] > int(args["--chunk"]) + int(args["--overlap"])
        ):
            chunk_by_busco(seq, seqs, busco_windows, args)
        else:
            if "--bed" in args and args["--bed"] and args["--bed"] != "None":
                stats = seq_stats(seq["seq"])
                bed_data[seq["title"]].append(
                    {
                        "start": seq["start"],
                        "end": seq["start"] + len(seq["seq"]),
                        "stats": stats,
                    }
                )
            if check_for_unmasked_bases(seq["seq"]):
                seqs.append((seq))
    if bed_data:
        write_bedfiles(bed_data, args)
    if args["--out"] is not None and args["--out"] != "None":
        chunked = ""
        for seq in seqs:
            if seq["length"] >= int(args["--min-length"]):
                chunked += ">%s_-_%d\n" % (seq["title"], seq["start"])
                chunked += "%s\n" % seq["seq"]
        outfile = args["--out"]
        if outfile.startswith("."):
            outfile = f'{args["--in"]}{args["--out"]}'
        with open(outfile, "w") as ofh:
            ofh.writelines(chunked)


if __name__ == "__main__":
    main()
