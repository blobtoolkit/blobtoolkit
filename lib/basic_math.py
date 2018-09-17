#!/usr/bin/env python3
"""Test module to learn how to use continuous integration."""


def square(x):
    """
    Return the square of number.

    >>> square(2)
    4
    >>> square(-2)
    4
    """
    return x * x


if __name__ == '__main__':
    import doctest
    doctest.testmod()
