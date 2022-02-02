#!/usr/bin/env python3
"""Functions to run external programs."""

import gzip
import shlex
from pathlib import Path
from subprocess import PIPE, run


def seqtk_subseq(infile, ids, outfile):
    """Run seqtk subseq."""
    cmd = "seqtk subseq %s -" % infile
    file_open = open
    if Path(infile).suffix == '.gz':
        file_open = gzip.open
    with file_open(outfile, 'wt') as ofh:
        process = run(shlex.split(cmd),
                      stdout=PIPE,
                      stderr=PIPE,
                      input=ids,
                      encoding='ascii')
        ofh.write(process.stdout)
    return process.returncode


if __name__ == '__main__':
    import doctest
    doctest.testmod()
