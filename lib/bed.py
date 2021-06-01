#!/usr/bin/env python3
# pylint: skip-file

"""Parse BED files into Fields."""


import math
import sys
from collections import Counter, OrderedDict, defaultdict
from glob import glob
from itertools import zip_longest
from os import path
from pathlib import Path

from tolkein import tofile
from tqdm import tqdm

import file_io
from field import Identifier, MultiArray, Variable
from run_external import seqtk_subseq


def field_settings():
    """Field settings for bed file imports."""
    return {
        "gc": {
            "meta": {
                "preload": True,
                "name": "GC",
                "scale": "scaleLinear",
                "type": "variable",
                "datatype": "float",
            },
            "parents": [
                {
                    "id": "gc_stats",
                    "scale": "scaleLinear",
                    "type": "variable",
                    "datatype": "float",
                },
                "children",
            ],
            "plot_axis": "x",
        },
        "length": {
            "meta": {
                "preload": True,
                "name": "Length",
                "clamp": 1,
                "scale": "scaleLog",
                "datatype": "integer",
                "type": "variable",
            },
            "parents": [
                {
                    "id": "length_stats",
                    "scale": "scaleLog",
                    "datatype": "integer",
                    "type": "variable",
                },
                "children",
            ],
            "plot_axis": "z",
        },
        "position": {
            "meta": {
                "name": "Position",
                "scale": "scaleLinear",
                "datatype": "integer",
                "type": "variable",
                "range": [0, 1],
            },
            "parents": [
                {
                    "id": "position_stats",
                    "scale": "scaleLinear",
                    "datatype": "integer",
                    "type": "variable",
                },
                "children",
            ],
        },
        "proportion": {
            "meta": {
                "name": "Proportion",
                "scale": "scaleLinear",
                "datatype": "float",
                "type": "variable",
                "range": [0, 1],
            },
            "parents": [
                {
                    "id": "proportion_stats",
                    "scale": "scaleLinear",
                    "datatype": "float",
                    "type": "variable",
                },
                "children",
            ],
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
                "name": "N",
                "scale": "scaleLinear",
                "datatype": "float",
                "type": "variable",
            },
            "parents": [
                {
                    "id": "n_stats",
                    "scale": "scaleLinear",
                    "datatype": "float",
                    "type": "variable",
                },
                "children",
            ],
        },
        "masked": {
            "meta": {
                "name": "Masked",
                "scale": "scaleLinear",
                "datatype": "float",
                "type": "variable",
            },
            "parents": [
                {
                    "id": "masked_stats",
                    "scale": "scaleLinear",
                    "datatype": "float",
                    "type": "variable",
                },
                "children",
            ],
        },
        "ncount": {
            "meta": {
                "name": "N count",
                "scale": "scaleLinear",
                "datatype": "integer",
                "type": "variable",
            },
            "parents": [
                {
                    "id": "ncount_stats",
                    "scale": "scaleLinear",
                    "datatype": "integer",
                    "type": "variable",
                },
                "children",
            ],
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


def parse_window_bed(filename):
    """Parse a BED file containing multiple values per sequence."""
    parsed = defaultdict(list)
    with tofile.open_file_handle(filename) as fh:
        for line in fh:
            row = line.strip().split("\t")
            parsed[row[0]].append((int(row[1]), float(row[4])))
    windowed = {}
    for seq_id, arr in parsed.items():
        windowed[seq_id] = [tup[1] for tup in sorted(arr, key=lambda tup: tup[0])]
    return windowed


def parse_full_tsv(filename):
    """Parse a TSV file containing one value per sequence."""
    values = defaultdict(dict)
    sd = defaultdict(dict)
    header = None
    with tofile.open_file_handle(filename) as fh:
        for line in fh:
            row = line.strip().split("\t")
            if header is None:
                header = {key: idx + 3 for idx, key in enumerate(row[3:])}
                continue
            values["length"][row[0]] = int(row[2])
            for key, idx in header.items():
                if key.endswith("_sd"):
                    sd[key.replace("_sd", "")][row[0]] = float(row[idx])
                else:
                    values[key][row[0]] = float(row[idx])
    return values, sd


def parse_windowed_tsv(filename):
    """Parse a TSV file containing one value per sequence."""
    values = defaultdict(lambda: defaultdict(list))
    sd = defaultdict(lambda: defaultdict(list))
    header = None
    with tofile.open_file_handle(filename) as fh:
        for line in fh:
            row = line.strip().split("\t")
            if header is None:
                header = {key: idx + 3 for idx, key in enumerate(row[3:])}
                continue
            values["length"][row[0]].append(int(row[2]))
            for key, idx in header.items():
                if key.endswith("_sd"):
                    sd[key.replace("_sd", "")][row[0]].append(float(row[idx]))
                else:
                    values[key][row[0]].append(float(row[idx]))
    return values, sd


def parse_tsvfiles(files):
    """Parse all tsvfiles."""
    windows = {}
    windows_sd = {}
    [fullfile, *files] = sorted(files, key=len)
    full, full_sd = parse_full_tsv(fullfile)
    for file in files:
        window = "".join(y for x, y in zip_longest(fullfile, file) if x != y).replace(
            ".tsv", ""
        )
        windows[window], windows_sd[window] = parse_windowed_tsv(file)
    return fullfile, windows, full


def parse_bedfiles(files):
    """Parse all bedfiles."""
    full = {}
    windows = {"0.1": {}}
    filenames = {}
    for filename in files:
        filepath = Path(filename)
        field = filepath.name.split(".")[-2]
        filenames[field] = filepath.name
        if field.endswith("_windows"):
            windows["0.1"][field.replace("_windows", "")] = parse_window_bed(filename)
        else:
            full[field] = parse_full_bed(filename)
    return filenames, windows, full


def parse(files, **kwargs):
    if "--bedtsvdir" in kwargs or "--bedtsvdir" in kwargs:
        if isinstance(files, str) and path.isdir(files):
            print("Reading all TSV files in %s" % files)
            files = glob("%s/*.tsv" % files)
        filename, all_windows, full = parse_tsvfiles(files)
        filenames = {"all": filename}
    else:
        if isinstance(files, str) and path.isdir(files):
            print("Reading all BED files in %s" % files)
            files = glob("%s/*.bed" % files)
        filenames, all_windows, full = parse_bedfiles(files)
    parsed = []
    settings = field_settings()
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
        key: {"range": [math.inf, -math.inf], "meta": {}} for key in settings.keys()
    }
    for field, data in full.items():
        filename = filenames.get(field, filenames.get("all", ""))
        if data:
            values = []
            for seq_id in identifiers.values:
                values.append(data[seq_id] if seq_id in data else 0)
            if values:
                meta = {}
                parents = []
                suffix = field.split("_")[-1]
                if suffix in settings:
                    meta = {
                        **settings[suffix]["meta"],
                        "field_id": field,
                        "file": filename,
                    }
                    if meta["datatype"] == "integer":
                        values = [int(value) for value in values]
                    value_range = [min(values), max(values)]
                    if "clamp" in meta and value_range[0] >= meta["clamp"]:
                        meta["clamp"] = False
                    parent_range = False
                    if "parents" in settings[suffix]:
                        parents = settings[suffix]["parents"]
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
                        if "range" in meta:
                            meta["range"][0] = min(meta["range"][0], value_range[0])
                            meta["range"][1] = max(meta["range"][1], value_range[1])
                        else:
                            meta["range"] = value_range
                    if "preload" in meta and meta["preload"] == 1:
                        if value_range[1] > ranges[suffix]["range"][1]:
                            meta["preload"] = True
                            if "plot_axis" in settings[suffix]:
                                kwargs["meta"].plot.update(
                                    {settings[suffix]["plot_axis"]: field}
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
                for window, windows in all_windows.items():
                    if field in windows:
                        window_values = []
                        for seq_id in identifiers.values:
                            values = windows[field][seq_id] if seq_id in data else []
                            if meta["datatype"] == "integer":
                                values = [[int(value)] for value in values]
                            else:
                                values = [[value] for value in values]
                            window_values.append(values)
                        parsed.append(
                            MultiArray(
                                "%s_windows_%s" % (field, window),
                                meta={
                                    "field_id": "%s_windows_%s" % (field, window),
                                    "name": "%s windows %s" % (meta["name"], window),
                                    "type": "multiarray",
                                    "datatype": meta["datatype"],
                                },
                                values=window_values,
                                parents=parents,
                                headers=[field],
                            )
                        )
    return parsed


def parent():
    """Set standard metadata."""
    return []
