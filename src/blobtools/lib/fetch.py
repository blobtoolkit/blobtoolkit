#!/usr/bin/env python3
"""Functions to fetch data and metadata."""

import glob
import os
from pathlib import Path

from ..lib import file_io
from .dataset import Metadata
from .field import Array
from .field import Category
from .field import Field
from .field import Identifier
from .field import MultiArray
from .field import Variable
from .taxdump import Taxdump

TYPES = {
    "identifier": Identifier,
    "category": Category,
    "variable": Variable,
    "array": Array,
    "multiarray": MultiArray,
    "field": Field,
}


def fetch_field(path_to_dataset, field_id, meta=None):
    """
    Load fields from file.

    fetch_field('tests/files/dataset', 'identifiers', meta)
    """
    field_meta = meta.field_meta(field_id)
    try:
        data = file_io.load_yaml("%s/%s.json" % (path_to_dataset, field_id))
        if data is None:
            data = file_io.load_yaml("%s/%s.json.gz" % (path_to_dataset, field_id))
        if data is not None:
            data.update({"meta": field_meta})
        field = TYPES[field_meta["type"]](field_id, **data)
    except TypeError:
        field = False
    except KeyError:
        field = False
    return field


def fetch_metadata(path_to_dataset, **kwargs):
    """
    Load Metadata from file.

    fetch_metadata('tests/files/dataset')
    """
    dataset_id = path_to_dataset.split("/").pop()
    new_meta = {}
    meta = None
    if not os.path.exists(path_to_dataset):
        os.makedirs(path_to_dataset)
    if kwargs.get("--meta"):
        new_meta = file_io.load_yaml(kwargs["--meta"])
        if (kwargs["--bed"] or kwargs["--fasta"]) and kwargs["--replace"]:
            files = glob.glob("%s/*" % kwargs["DIRECTORY"])
            for file in files:
                os.remove(file)
    try:
        meta = kwargs["meta"]
    except KeyError:
        try:
            meta = file_io.load_yaml("%s/meta.json" % path_to_dataset)
            if meta is None:
                meta = file_io.load_yaml("%s/meta.json.gz" % path_to_dataset)
        except ValueError:
            pass
    if meta is None:
        meta = {}
    if "id" not in meta:
        meta["id"] = dataset_id
        meta["name"] = dataset_id
    for key, value in new_meta.items():
        if isinstance(value, dict):
            try:
                meta[key].update({k: v for k, v in value.items()})
            except KeyError:
                meta[key] = value
        elif isinstance(value, list):
            meta[key] += value
        else:
            meta[key] = value
    return Metadata(dataset_id, **meta)


def fetch_taxdump(path_to_taxdump):
    """Load Taxdump from file."""
    json_file = "%s/taxdump.json" % path_to_taxdump
    if not Path(json_file).exists():
        print("Parsing taxdump")
    else:
        print("Loading parsed taxdump")
    data = file_io.load_yaml(json_file)
    if data is None:
        taxdump = Taxdump(path_to_taxdump)
        file_io.write_file(json_file, taxdump.values_to_dict())
    else:
        taxdump = Taxdump(path_to_taxdump, **data)
    return taxdump


if __name__ == "__main__":
    import doctest

    doctest.testmod()
