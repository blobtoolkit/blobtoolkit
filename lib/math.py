#!/usr/bin/env python3
"""
Example module with basic maths functions.

Tests are included in docstrings and run using doctest.
"""


def square(x):
    """
    Return the square of number.

    >>> square(2)
    4
    >>> square(-2)
    4
    """
    return x * x


def cube(x):
    """
    Return the cube of number.

    >>> cube(2)
    8
    >>> cube(-2)
    -8
    """
    return x ** 3


if __name__ == '__main__':
    import doctest
    doctest.testmod()
