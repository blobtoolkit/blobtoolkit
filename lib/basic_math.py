#!/usr/bin/env python3
"""Test module to learn how to use continuous integration."""


def square(number):
    """
    Return the square of number.

    >>> square(2)
    4
    >>> square(-2)
    4
    """
    return number * number


if __name__ == '__main__':
    import doctest
    doctest.testmod()
