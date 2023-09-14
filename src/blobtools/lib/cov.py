#!/usr/bin/env python3

# pylint: disable=no-member, too-many-locals
# pylint: skip-file

"""Parse base and read coverage into Variable Field."""

import math
import os
import sys
from collections import defaultdict
from pathlib import Path

from blobtk import depth
from blobtk import filter
from tqdm import tqdm

from .field import Variable
from .file_io import load_yaml


def get_coverage(bam_file):
    """Get base coverage."""
    binned_covs = depth.bam_to_depth(bam=bam_file)
    return {cov.seq_name: cov.bins[0] for cov in binned_covs}


def parse_bam(bam_file, **kwargs):
    """Parse coverage into a Variables."""
    identifiers = kwargs["dependencies"]["identifiers"]
    ids = identifiers.values
    lengths = kwargs["dependencies"]["length"].values
    ncounts = kwargs["dependencies"]["ncount"].values
    parts = bam_file.split("=")
    base_name = parts[1]
    f_char = Path(parts[0]).suffix[1]
    index_file = Path(f"{parts[0]}.csi")
    if index_file.is_file():
        index_file = False
    print(f"Loading mapping data from {parts[0]} as {parts[1]}")
    _covs = get_coverage(parts[0])
    if index_file:
        os.remove(index_file)
    if not identifiers.validate_list(list(_covs.keys())):
        raise UserWarning(
            "Contig names in the coverage file did not match dataset identifiers."
        )
    covs = [float("%.4f" % (_covs[seq_id]) if seq_id in _covs else 0) for seq_id in ids]
    field_id = f"{base_name}_cov"
    fields = {
        "cov_id": field_id,
        "cov_range": [
            min(covs + [kwargs["cov_range"][0]]),
            max(covs + [kwargs["cov_range"][1]]),
        ],
    }
    fields["cov"] = Variable(
        field_id,
        values=covs,
        meta={"field_id": field_id, "file": bam_file},
        parents=[
            "children",
            {
                "id": "base_coverage",
                "clamp": 0.01 if fields["cov_range"][0] == 0 else False,
                "range": fields["cov_range"],
            },
            "children",
        ],
    )
    return fields


def apply_filter(ids, fastq_files, **kwargs):
    """Filter FASTQ file based on read alignment file."""
    suffix = kwargs["--suffix"]
    bam_file = kwargs["--cov"]
    options = {"list": ids, "suffix": suffix, "fastq_out": True}
    filetype = Path(bam_file).suffix
    if filetype == ".bam":
        options["bam"] = bam_file
    elif filetype == ".cram":
        options["cram"] = bam_file
    else:
        print("ERROR: Alignment file suffix %s is not suported", filetype)
        sys.exit(1)
    options["fastq1"] = fastq_files[0]
    if len(fastq_files) > 1:
        options["fastq2"] = fastq_files[1]
    filter.fastx(**options)


def parse_json_cov(json_file, **kwargs):
    """Parse coverage from JSON cov file."""
    parts = json_file.split("=")
    base_name = parts[1]
    data = load_yaml(parts[0])
    covs = []
    if "values" in data:
        for value in data["values"]:
            covs.append(float("%.4f" % value))
    if base_name.endswith("_read_cov"):
        type = "read_cov"
        parent = "read_coverage"
        datatype = "float"
        clamp = 1
    elif base_name.endswith("_cov"):
        type = "cov"
        parent = "base_coverage"
        datatype = "integer"
        clamp = 0.01
    else:
        return None
    field_id = base_name
    fields = {}
    fields["%s_id" % type] = field_id
    fields["%s_range" % type] = [
        min(covs + [kwargs["%s_range" % type][0]]),
        max(covs + [kwargs["%s_range" % type][1]]),
    ]
    if kwargs["meta"].has_field(field_id):
        file_name = kwargs["meta"].field_meta(field_id)["file"]
    else:
        file_name = json_file
    fields[type] = Variable(
        field_id,
        values=covs,
        meta={"field_id": field_id, "file": file_name},
        parents=[
            "children",
            {
                "id": parent,
                "datatype": "integer",
                "clamp": clamp if fields["%s_range" % type][0] == 0 else False,
                "range": fields["%s_range" % type],
            },
            "children",
        ],
    )
    return fields


def base_names(files):
    """Set file base names."""
    names = {}
    unique = True
    for file in files:
        if "=" in file:
            parts = file.split("=")
            name = parts[1]
            names[name] = parts[0]
        else:
            name = Path(file).stem
            if name in names:
                unique = False
            else:
                names[name] = file
    if not unique:
        print("ERROR: Unable to set unique names for coverage files")
        sys.exit(1)
    return ["%s=%s" % (v, k) for k, v in names.items()]


def parse(files, **kwargs):
    """Parse all BAM files."""
    parsed = []
    if kwargs["meta"].has_field("base_coverage"):
        cov_range = kwargs["meta"].field_meta("base_coverage")["range"]
    else:
        cov_range = [math.inf, -math.inf]
    if kwargs["meta"].has_field("read_coverage"):
        read_cov_range = kwargs["meta"].field_meta("read_coverage")["range"]
    else:
        read_cov_range = [math.inf, -math.inf]
    names = base_names(files)
    for file in names:
        if ".json" in file:
            fields = parse_json_cov(
                file, **kwargs, cov_range=cov_range, read_cov_range=read_cov_range
            )
        else:
            fields = parse_bam(
                file, **kwargs, cov_range=cov_range, read_cov_range=read_cov_range
            )
        if "cov" in fields:
            parsed.append(fields["cov"])
            cov_range = fields["cov_range"]
            if "y" not in kwargs["meta"].plot:
                kwargs["meta"].plot.update({"y": fields["cov_id"]})
        if "read_cov" in fields:
            parsed.append(fields["read_cov"])
            read_cov_range = fields["read_cov_range"]
    return parsed


def parent():
    """Set standard metadata for Coverage."""
    coverage = {
        "datatype": "float",
        "type": "variable",
        "scale": "scaleLog",
        "id": "coverage",
        "name": "coverage",
    }
    return [coverage]


def weighted_mean(values, weights, log=False):
    """Calculate weighted mean and standard deviation."""
    if log:
        mean = sum(
            [
                math.log10(value + 0.01) * weight
                for value, weight in zip(values, weights)
            ]
        ) / sum(weights)
        mean = 10**mean - 0.01
    else:
        mean = sum([value * weight for value, weight in zip(values, weights)]) / sum(
            weights
        )
    return mean


def summarise(indices, fields, **kwargs):
    """Summarise coverage."""
    summary = {}
    for library in fields["libraries"]:
        stats = {}
        covs = [fields[f"{library}_cov"].values[i] for i in indices]
        lengths = [fields["length"].values[i] for i in indices]
        coverage = weighted_mean(covs, lengths, log=True)
        if library in kwargs["meta"].reads:
            try:
                platform = kwargs["meta"].reads[library]["platform"]
            except KeyError:
                platform = "unknown"
            strategy = "single"
            try:
                if "strategy" in kwargs["meta"].reads[library]:
                    strategy = kwargs["meta"].reads[library]["strategy"]
                elif (
                    "file" in kwargs["meta"].reads[library]
                    and ";" in kwargs["meta"].reads[library]["file"]
                ):
                    strategy = "paired"
            except KeyError:
                continue
            stats.update(
                {
                    "coverage": float("%.3f" % coverage),
                    "platform": platform,
                    "strategy": strategy,
                }
            )
            if "file" in kwargs["meta"].reads[library]:
                stats["file"] = kwargs["meta"].reads[library]["file"].split(";")
            if "url" in kwargs["meta"].reads[library]:
                stats["url"] = kwargs["meta"].reads[library]["url"]
        else:
            stats["coverage"] = float("%.3f" % coverage)
        summary[library] = stats
    return summary


def remove_from_meta(meta):
    """Delete all coverage fields."""
    field_ids = []
    if meta.has_field("coverage"):
        field_ids += meta.remove_field("coverage")
        meta.plot.pop("y", None)
    meta.reads = {}
    return field_ids
