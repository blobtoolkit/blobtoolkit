#!/usr/bin/env python3
"""
Example module with functions to read, write and parse files.

Simple tests are included in docstrings and run using doctest.

Further tests are implemented in a separate file using pytest
(tests/unit_tests/test_file_io.py) to avoid adding test-specific code to
docstrings.
"""

import gzip
import json
import pathlib
import yaml


def read_file(filename):
    r"""
    Read a file, automatically detect gzipped filed based on suffix.

    >>> read_file('tests/files/infile')
    'testfile content\n'
    >>> read_file('nofile')
    ''
    """
    if '.gz' in pathlib.Path(filename).suffixes:
        try:
            with gzip.open(filename, 'rt') as fh:
                return fh.read()
        except OSError:
            return ''
    try:
        with open(filename, 'r') as fh:
            return fh.read()
    except IOError:
        return ''


def load_yaml(filename):
    """
    Parse a JSON/YAML file.

    load_yaml('identifiers.yaml')
    """
    data = read_file(filename)
    identifiers = yaml.load(data)
    return identifiers


def write_file(filename, data):
    """
    Write a file, use suffix to determine type and compression.

    - types: '.json', '.yaml'
    - compression: None, '.gz'

    write_file('variable.json.gz')
    """
    if '.json' in filename:
        content = json.dumps(data, indent=1)
    elif '.yaml' in filename:
        content = yaml.dump(data, indent=1)
    else:
        content = data
    if '.gz' in filename:
        try:
            with gzip.open(filename, 'wt') as fh:
                fh.write(content)
        except OSError:
            return False
    else:
        try:
            with open(filename, 'wt') as fh:
                fh.write(content)
        except IOError:
            return False
    return True


if __name__ == '__main__':
    import doctest
    doctest.testmod()
