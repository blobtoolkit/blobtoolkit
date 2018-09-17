#!/usr/bin/env python3
"""functions to read, write and parse files."""


def read_file(filename):
    """
    Reads a file.

    >>> read_file('test/files/infile')
    'testfile content\\n'
    """
    with open(filename, 'r') as fh:
        return fh.read()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
