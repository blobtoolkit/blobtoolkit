#!/usr/bin/env python3
# pylint: skip-file

"""Parse BED files into Fields."""


import math
import sys
from collections import Counter
from collections import OrderedDict
from collections import defaultdict
from copy import deepcopy
from glob import glob
from itertools import zip_longest
from os import path
from os import stat
from pathlib import Path

from tolkein import tofile
from tqdm import tqdm

from .field import Array
from .field import Identifier
from .field import MultiArray
from .field import Variable


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
                    "id": "gc_data",
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
                    "id": "length_data",
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
                    "id": "position_data",
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
                    "id": "proportion_data",
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
                "range": [math.inf, -math.inf],
            },
            "parents": [
                {
                    "id": "coverage",
                    "type": "variable",
                    "datatype": "float",
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
                    "id": "n_data",
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
                    "id": "masked_data",
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
                    "id": "ncount_data",
                    "scale": "scaleLinear",
                    "datatype": "integer",
                    "type": "variable",
                },
                "children",
            ],
        },
        "count": {
            "meta": {
                "scale": "scaleLinear",
                "name": "Count",
                "datatype": "integer",
                "type": "variable",
            },
            "parents": [
                {
                    "id": "count_data",
                    "type": "variable",
                    "datatype": "integer",
                },
                "children",
            ],
        },
        "cpm": {
            "meta": {
                "scale": "scaleLinear",
                "name": "Count per Mb",
                "datatype": "float",
                "type": "variable",
            },
            "parents": [
                {
                    "id": "cpm_data",
                    "type": "variable",
                    "datatype": "float",
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
    n = defaultdict(dict)
    header = None
    with tofile.open_file_handle(filename) as fh:
        for line in fh:
            row = line.strip().split("\t")
            if header is None:
                header = {key: idx + 3 for idx, key in enumerate(row[3:])}
                continue
            length = int(row[2]) - int(row[1])
            values["length"][row[0]] = length
            values["position"][row[0]] = int(row[2])
            for key, idx in header.items():
                if key.endswith("_sd"):
                    sd[key[: key.rfind("_sd")]][row[0]] = float(row[idx])
                elif key.endswith("_n"):
                    n[key[: key.rfind("_n")]][row[0]] = float(row[idx])
                elif key.endswith("_cpm"):
                    values[key][row[0]] = float(
                        "%.3g" % (float(row[idx]) / length * 1000000)
                    )
                elif key.endswith("_count"):
                    values[key][row[0]] = int(row[idx])
                else:
                    values[key][row[0]] = float(row[idx])
    return values, sd, n


def parse_windowed_tsv(filename, window):
    """Parse a TSV file containing one value per sequence."""
    values = defaultdict(lambda: defaultdict(list))
    sd = defaultdict(lambda: defaultdict(list))
    n = defaultdict(lambda: defaultdict(list))
    lengths = {}
    header = None
    with tofile.open_file_handle(filename) as fh:
        for line in fh:
            row = line.strip().split("\t")
            if header is None:
                header = {key: idx + 3 for idx, key in enumerate(row[3:])}
                continue
            length = int(row[2]) - int(row[1])
            if float(window) > 1:
                prev_len = lengths.get(row[0], 0)
                if length < prev_len:
                    continue
                lengths[row[0]] = length
            values["length"][row[0]].append(length)
            values["position"][row[0]].append(round(int(row[1]) + length / 2))
            for key, idx in header.items():
                if key.endswith("_sd"):
                    sd[key[: key.rfind("_sd")]][row[0]].append(float(row[idx]))
                elif key.endswith("_n"):
                    n[key[: key.rfind("_n")]][row[0]].append(float(row[idx]))
                elif key.endswith("_cpm"):
                    values[key][row[0]].append(
                        float("%.3g" % (float(row[idx]) / length * 1000000))
                    )
                elif key.endswith("_count"):
                    values[key][row[0]].append(int(row[idx]))
                else:
                    values[key][row[0]].append(float(row[idx]))
    return values, sd, n


def parse_tsvfiles(files):
    """Parse all tsvfiles."""
    windows = {}
    windows_sd = {}
    windows_n = {}
    [fullfile, *files] = sorted(files, key=len)
    full, sd, n = parse_full_tsv(fullfile)
    for file in files:
        if stat(file).st_size == 0:
            continue
        window = "".join(y for x, y in zip_longest(fullfile, file) if x != y).replace(
            ".tsv", ""
        )
        windows[window], windows_sd[window], windows_n[window] = parse_windowed_tsv(
            file, window
        )
    return (
        fullfile,
        {"values": windows, "sd": windows_sd, "n": windows_n},
        {"values": full, "sd": sd, "n": n},
    )


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
            windows["0.1"][field[: field.rfind("_windows")]] = parse_window_bed(
                filename
            )
        else:
            full[field] = parse_full_bed(filename)
    return filenames, windows, full


def validate_range(meta):
    """Ensure range is greater than zero."""
    if meta["range"][1] < meta["range"][0]:
        return False
    if meta["range"][1] == meta["range"][0]:
        if meta["datatype"] == "integer":
            meta["range"][1] += 1
        elif meta["range"][1] == 0:
            meta["range"][1] += 0.1
        else:
            mag = math.floor(math.log10(meta["range"][1]))
            meta["range"][1] += 10 ** (mag - 1)
    return True


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
    full_n = {}
    full_sd = {}
    if isinstance(full, dict):
        full_sd = full["sd"]
        full_n = full["n"]
        full = full["values"]
    all_windows_n = {}
    all_windows_sd = {}
    if isinstance(all_windows, dict):
        all_windows_n = all_windows["n"]
        all_windows_sd = all_windows["sd"]
        all_windows = all_windows["values"]
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
                    meta = deepcopy(settings[suffix]["meta"])
                    meta.update(
                        {
                            "field_id": field,
                            "file": filename,
                        }
                    )
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
                    if not validate_range(meta):
                        continue
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
                    if field.endswith("_%s" % suffix):
                        meta["name"] = "%s %s" % (
                            field[: field.rfind("_%s" % suffix)],
                            meta["name"],
                        )
                if not meta:
                    continue
                parsed.append(
                    Variable(
                        field,
                        meta=meta,
                        values=values,
                        parents=parents,
                    )
                )
                if field in full_sd:
                    stats_values = []
                    for seq_id in identifiers.values:
                        values = (
                            [full_sd[field][seq_id], full_n[field][seq_id]]
                            if seq_id in data
                            else []
                        )
                        stats_values.append(values)
                    parsed.append(
                        Array(
                            "%s_stats" % field,
                            meta={
                                "field_id": "%s_stats" % field,
                                "name": "%s stats" % meta["name"],
                                "type": "array",
                                "datatype": "mixed",
                            },
                            values=stats_values,
                            parents=parents,
                            headers=["sd", "n"],
                        )
                    )
                for window, windows in all_windows.items():
                    windows_sd = all_windows_sd.get(window, {})
                    windows_n = all_windows_n.get(window, {})
                    if field in windows:
                        window_values = []
                        headers = [field]
                        if field in windows_sd:
                            headers += ["sd", "n"]
                        for seq_id in identifiers.values:
                            seq_values = []
                            if seq_id in data:
                                for idx, value in enumerate(windows[field][seq_id]):
                                    if meta["datatype"] == "integer":
                                        value = int(value)
                                    if field in windows_sd:
                                        value = [
                                            value,
                                            windows_sd[field][seq_id][idx],
                                            windows_n[field][seq_id][idx],
                                        ]
                                    else:
                                        value = [value]
                                    seq_values.append(value)
                            window_values.append(seq_values)
                        windows_field = "%s_windows" % field
                        if str(window) != "0.1":
                            windows_field += "_%s" % str(window)
                        parsed.append(
                            MultiArray(
                                windows_field,
                                meta={
                                    "field_id": windows_field,
                                    "name": "%s windows %s" % (meta["name"], window),
                                    "type": "multiarray",
                                    "datatype": "mixed",
                                },
                                values=window_values,
                                parents=parents,
                                headers=headers,
                            )
                        )
    return parsed


def parent():
    """Set standard metadata."""
    return []
