#!/usr/bin/env python3

# pylint: disable=c-extension-no-member

"""
Read, write and parse files.

Simple tests are included in docstrings and run using doctest.

Further tests are implemented in a separate file using pytest
(tests/unit_tests/test_file_io.py) to avoid adding test-specific code to
docstrings.
"""

import csv
import gzip
import io
import os
import pathlib
import sys
from itertools import groupby
from subprocess import PIPE
from subprocess import Popen

import ujson
import yaml


def read_file(filename):
    r"""
    Read a whole file into memory.

    >>> read_file('tests/files/infile')
    'testfile content\n'
    >>> read_file('nofile')

    """
    try:
        with stream_file(filename) as fh:
            return fh.read()
    except AttributeError:
        return None


def stream_file(filename):
    """
    Stream a file, line by line.

    Automatically detect gzipped files based on suffix.
    """
    if ".gz" in pathlib.Path(filename).suffixes:
        try:
            return gzip.open(filename, "rt")
        except OSError:
            return None
    try:
        return open(filename, "r")
    except IOError:
        return None


def stream_fasta(filename):
    """
    Stream a fasta file, sequence by sequence.

    Automatically detect gzipped files based on suffix.
    """
    if ".gz" in pathlib.Path(filename).suffixes:
        cmd = ["gunzip", "-c", filename]
    else:
        cmd = ["cat", filename]
    with Popen(cmd, stdout=PIPE, encoding="utf-8", bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == ">"))
        for header in faiter:
            seq_id = header.__next__().split()[0].replace(">", "")
            seq_str = "".join(map(lambda s: s.strip(), faiter.__next__()))
            yield seq_id, seq_str


def load_yaml(filename):
    """
    Parse a JSON/YAML file.

    load_yaml('identifiers.yaml')
    """
    data = read_file(filename)
    if data is None:
        return data
    if ".json" in filename:
        content = ujson.loads(data)
    elif ".yaml" in filename:
        content = yaml.full_load(data)
    else:
        content = data
    return content


def write_file(filename, data, plain=False):  # pylint: disable=too-many-branches
    """
    Write a file, use suffix to determine type and compression.

    - types: '.json', '.yaml'
    - compression: None, '.gz'

    write_file('variable.json.gz')
    """
    if ".json" in filename:
        content = ujson.dumps(data, indent=1, escape_forward_slashes=False)
    elif ".yaml" in filename:
        content = yaml.dump(data, indent=1)
    elif filename == "STDOUT":
        sys.stdout.write(
            ujson.dumps(data, indent=1, escape_forward_slashes=False) + "\n"
        )
        return True
    elif filename == "STDOUT":
        sys.stderr.write(
            ujson.dumps(data, indent=1, escape_forward_slashes=False) + "\n"
        )
        return True
    elif plain:
        content = "\n".join(data)
    elif ".csv" in filename or ".tsv" in filename:
        output = io.StringIO()
        if ".csv" in filename:
            writer = csv.writer(output, quoting=csv.QUOTE_NONNUMERIC)
        else:
            writer = csv.writer(output, delimiter="\t")
        for row in data:
            writer.writerow(row)
        content = output.getvalue()
    else:
        content = data
    if ".gz" in filename:
        try:
            with gzip.open(filename, "wt") as fh:
                fh.write(content)
        except OSError:
            return False
    else:
        try:
            with open(filename, "wt") as fh:
                fh.write(content)
        except IOError:
            return False
    return True


def delete_file(file):
    """Delete a file if exists."""
    if os.path.exists(file):
        os.remove(file)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
