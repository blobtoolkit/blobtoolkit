#!/usr/bin/env python3

"""
Generate plots using BlobToolKit Viewer.
"""

from blobtoolkit_view import view


def cli():
    view.cli("blobtools view")


if __name__ == "__main__":
    cli()
