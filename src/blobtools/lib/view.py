#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-locals, too-many-statements, broad-except

"""
Generate plots using BlobToolKit Viewer.
"""

from blobtoolkit_view import view


def cli():
    view.cli("blobtools view")


if __name__ == "__main__":
    cli()
