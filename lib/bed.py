#!/usr/bin/env python3
# pylint: skip-file

"""Parse BED files into Fields."""

import math
import sys
from collections import Counter, OrderedDict
from glob import glob
from os import path
from pathlib import Path

from tolkein import tofile
from tqdm import tqdm

import file_io
from field import Identifier, Variable
from run_external import seqtk_subseq

SETTINGS = {
    "gc": {
        "meta": {
            "preload": True,
            "scale": "scaleLinear",
            "name": "GC",
            "datatype": "float",
            "type": "variable",
        }
    },
    "length": {
        "meta": {
            "preload": True,
            "scale": "scaleLog",
            "name": "Length",
            "clamp": 1,
            "datatype": "integer",
            "type": "variable",
        }
    },
    "cov": {
        "meta": {
            "preload": 1,
            "scale": "scaleLog",
            "name": "Coverage",
            "clamp": 0.01,
            "datatype": "float",
            "type": "variable",
        },
        "parents": [
            {
                "id": "coverage",
                "type": "variable",
                "datatype": "float",
                "range": [math.inf, -math.inf],
            },
            "children",
        ],
        "plot_axis": "y",
    },
    "n": {
        "meta": {
            "scale": "scaleLinear",
            "name": "N",
            "datatype": "float",
            "type": "variable",
        }
    },
    "masked": {
        "meta": {
            "scale": "scaleLinear",
            "name": "Masked",
            "datatype": "float",
            "type": "variable",
        }
    },
    "ncount": {
        "meta": {
            "scale": "scaleLinear",
            "name": "N count",
            "datatype": "integer",
            "type": "variable",
        }
    },
}


def parse_full_bed(filename):
    """Parse a BED file containing one value per sequence."""
    parsed = {}
    with tofile.open_file_handle(filename) as fh:
        for line in fh:
            row = line.strip().split("\t")
            if len(row) < 5:
                return parsed
            parsed[row[0]] = float(row[4])
    return parsed


def parse(files, **kwargs):
    if isinstance(files, str) and path.isdir(files):
        print("Reading all BED files in %s" % files)
        files = glob("%s/*.bed" % files)
    parsed = []
    full = {}
    windows = {}
    filenames = {}
    for filename in files:
        filepath = Path(filename)
        field = filepath.name.split(".")[-2]
        filenames[field] = filepath.name
        if field.endswith("_windows"):
            windows[field.replace("_windows", "")] = filename
        else:
            full[field] = parse_full_bed(filename)
    identifiers = kwargs["dependencies"]["identifiers"]
    keys = []
    if "length" in full:
        keys = list(full["length"].keys())
        lengths = list(full["length"].values())
        kwargs["meta"].assembly.update({"span": sum(lengths)})
        if "z" not in kwargs["meta"].plot:
            kwargs["meta"].plot.update({"z": "length"})
    if "gc" in full and "x" not in kwargs["meta"].plot:
        kwargs["meta"].plot.update({"x": "gc"})
    if not identifiers:
        if not keys:
            print("ERROR: Unable to set identifiers")
            sys.exit(1)
        identifiers = Identifier(
            "identifiers",
            meta={"field_id": "identifiers"},
            values=keys,
            parents=[],
        )
        kwargs["meta"].assembly.update({"scaffold-count": len(identifiers.values)})
        parsed.append(identifiers)
    ranges = {
        key: {"range": [math.inf, -math.inf], "meta": {}} for key in SETTINGS.keys()
    }
    for field, data in full.items():
        if data:
            values = []
            for seq_id in identifiers.values:
                values.append(data[seq_id] if seq_id in data else 0)
            if values:
                meta = {}
                parents = []
                suffix = field.split("_")[-1]
                if suffix in SETTINGS:
                    meta = {
                        **SETTINGS[suffix]["meta"],
                        "field_id": field,
                        "file": filenames[field],
                    }
                    if meta["datatype"] == "integer":
                        values = [int(value) for value in values]
                    value_range = [min(values), max(values)]
                    if "clamp" in meta and value_range[0] >= meta["clamp"]:
                        meta["clamp"] = False
                    parent_range = False
                    if "parents" in SETTINGS[suffix]:
                        parents = SETTINGS[suffix]["parents"]
                        for parent in parents:
                            if "range" in parent:
                                parent_range = True
                                parent["range"][0] = min(
                                    parent["range"][0], value_range[0]
                                )
                                parent["range"][1] = max(
                                    parent["range"][1], value_range[1]
                                )
                    if not parent_range:
                        meta["range"] = value_range
                    if "preload" in meta and meta["preload"] == 1:
                        if value_range[1] > ranges[suffix]["range"][1]:
                            meta["preload"] = True
                            if "plot_axis" in SETTINGS[suffix]:
                                kwargs["meta"].plot.update(
                                    {SETTINGS[suffix]["plot_axis"]: field}
                                )
                            if "preload" in ranges[suffix]["meta"]:
                                ranges[suffix]["meta"]["preload"] = False
                            ranges[suffix].update({"range": value_range, "meta": meta})
                        else:
                            meta["preload"] = False
                parsed.append(
                    Variable(
                        field,
                        meta=meta,
                        values=values,
                        parents=parents,
                    )
                )
    return parsed


def parent():
    """Set standard metadata."""
    return []
