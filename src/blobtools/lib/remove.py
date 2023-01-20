#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-nested-blocks

"""
Remove fields from a BlobDir.

Usage:
    blobtools remove [--all] [--busco] [--cov] [--fasta] [--field STRING...] [--hits]
                     DIRECTORY

Arguments:
    DIRECTORY             Existing Blob directory.

Options:
    --all             Remove all fields except identifiers.
    --busco           Remove all BUSCO fields.
    --cov             Remove all cov and read_cov fields.
    --fasta           Remove gc, length and ncount fields.
    --field STRING    Remove fields by ID.
    --hits            Remove all taxonomy fields.

Examples:
    # 1. Remove BUSCO and ncount fields from a BlobDir
    blobtools remove --busco --field ncount BlobDir

"""

import glob
import sys
from collections import defaultdict

from docopt import docopt

from ..lib import busco
from ..lib import cov
from ..lib import fasta
from ..lib import file_io
from ..lib import hits
from .fetch import fetch_metadata
from .version import __version__

FIELDS = [
    {"flag": "--fasta", "module": fasta},
    {"flag": "--busco", "module": busco},
    {"flag": "--cov", "module": cov},
    {"flag": "--hits", "module": hits},
]


def remove_field(meta, field_id):
    """Delete field and any descendants."""
    field_ids = []
    if meta.has_field(field_id):
        field_ids = meta.remove_field(field_id)
    meta.plot = {
        key: value for key, value in meta.plot.items() if value not in field_ids
    }
    return field_ids


def remove_read_metadata(meta, field_ids):
    """Delete any read metadata for remove fields."""
    covs = defaultdict(dict)
    for field_id in field_ids:
        if field_id.endswith("_cov"):
            if field_id.endswith("_read_cov"):
                root = field_id.replace("_read_cov", "")
                covs[root].update({"read": True})
            else:
                root = field_id.replace("_cov", "")
                covs[root].update({"base": True})
    for key, value in covs.items():
        if "base" in value and "read" in value:
            meta.reads.pop(key, None)


def remove_static_plots(meta, directory):
    """Delete static plots."""
    files = glob.glob("%s/*.??g" % directory)
    for file in files:
        file_io.delete_file(file)
    setattr(meta, "static_plots", False)


def main(args):
    """Entrypoint for blobtools remove."""
    meta = fetch_metadata(args["DIRECTORY"], **args)
    field_ids = []
    for field in FIELDS:
        if args[field["flag"]] or args["--all"]:
            field_ids += field["module"].remove_from_meta(meta=meta)
    for field_id in args["--field"]:
        field_ids += remove_field(meta, field_id)
    if args["--all"]:
        for field_id in meta.list_fields():
            if field_id != "identifiers":
                field_ids += remove_field(meta, field_id)
    if meta.reads:
        remove_read_metadata(meta, field_ids)
    if meta.plot:
        axes = ["x", "y", "z", "cat"]
        for axis in axes:
            if axis not in meta.plot:
                remove_static_plots(meta, args["DIRECTORY"])
                break
    else:
        remove_static_plots(meta, args["DIRECTORY"])
    if field_ids:
        file_io.delete_file("%s/CHECKSUM" % args["DIRECTORY"])
        file_io.delete_file("%s/summary.json" % args["DIRECTORY"])
        for field_id in field_ids:
            file_io.delete_file("%s/%s.json" % (args["DIRECTORY"], field_id))
        file_io.write_file("%s/meta.json" % args["DIRECTORY"], meta.to_dict())


def cli():
    """Entry point."""
    if len(sys.argv) == sys.argv.index(__name__.split(".")[-1]) + 1:
        args = docopt(__doc__, argv=[])
    else:
        args = docopt(__doc__, version=__version__)
    main(args)


if __name__ == "__main__":
    cli()
