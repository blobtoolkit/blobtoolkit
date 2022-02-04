#!/usr/bin/env python3
"""Convert a blobDB to BlobDir Fields."""

# pylint: disable=too-many-locals

import math
from collections import defaultdict
from pathlib import Path

from ..lib import cov
from ..lib import file_io
from ..lib import hits
from .field import Category
from .field import Identifier
from .field import MultiArray
from .field import Variable


def field_name_from_path(path):
    """Extract field name from file path."""
    parts = Path(path).stem.split(".")
    field_name = parts[-1]
    if len(parts) > 1:
        if parts[-1] in ("bam", "sam", "cram"):
            field_name = parts[-2]
    return field_name


def values_from_blob_db(blob_db):
    """Read values from a blobDB into a dict of lists of values."""
    values = defaultdict(list)
    for identifier in blob_db["order_of_blobs"]:
        blob = blob_db["dict_of_blobs"][identifier]
        values["lengths"].append(blob.get("length", 0))
        values["gcs"].append(blob.get("gc", 0))
        values["n_counts"].append(blob.get("n_count", 0))
        for cov_lib in blob_db["covLibs"].keys():
            values["%s_cov" % cov_lib].append(blob["covs"].get(cov_lib, 0))
            values["%s_read_cov" % cov_lib].append(blob["read_cov"].get(cov_lib, 0))
        for tax_rule in blob_db["taxrules"]:
            for rank, results in blob["taxonomy"][tax_rule].items():
                values["%s_%s" % (tax_rule, rank)].append(results.get("tax", "no-hit"))
                values["%s_%s_score" % (tax_rule, rank)].append(
                    float(results.get("score", 0))
                )
                values["%s_%s_cindex" % (tax_rule, rank)].append(
                    int(results.get("c_index", 0))
                )
    return values


def hits_from_blob_db(blob_db, tax_rule):
    """Read values from a blobDB into a dict of lists of values."""
    values = defaultdict(list)
    for identifier in blob_db["order_of_blobs"]:
        if identifier in values and tax_rule == "bestsumorder":
            continue
        hit_list = blob_db["dict_of_blobs"][identifier]["hits"]
        for tax_lib in blob_db["hitLibs"].keys():
            if tax_lib in hit_list:
                for detail in hit_list[tax_lib]:
                    try:
                        values[identifier].append(
                            {
                                "taxid": int(detail["taxId"]),
                                "score": int(detail["score"]),
                            }
                        )
                    except ValueError:
                        pass
    values = [values.get(identifier, []) for identifier in blob_db["order_of_blobs"]]
    return values


def parse(file, **kwargs):
    """Parse all synonym files."""
    blob_db = file_io.load_yaml(file)
    kwargs["meta"].assembly.update({"file": blob_db["assembly_f"]})
    parsed = []
    identifiers = kwargs["dependencies"]["identifiers"]
    if not identifiers:
        identifiers = Identifier(
            "identifiers",
            meta={"field_id": "identifiers"},
            values=blob_db["order_of_blobs"],
            parents=[],
        )
        kwargs["meta"].assembly.update({"scaffold-count": len(identifiers.values)})
        parsed.append(identifiers)
    values = values_from_blob_db(blob_db)
    kwargs["meta"].assembly.update({"span": sum(values["lengths"])})
    parsed.append(
        Variable(
            "gc",
            meta={
                "preload": True,
                "scale": "scaleLinear",
                "field_id": "gc",
                "name": "GC",
                "datatype": "float",
                "range": [min(values["gcs"]), max(values["gcs"])],
            },
            values=values["gcs"],
            parents=[],
        )
    )
    _min = min(values["lengths"])
    parsed.append(
        Variable(
            "length",
            meta={
                "field_id": "length",
                "preload": True,
                "scale": "scaleLog",
                "name": "Length",
                "clamp": 100 if _min == 0 else False,
                "datatype": "integer",
                "range": [_min, max(values["lengths"])],
            },
            parents=[],
            values=values["lengths"],
        )
    )
    parsed.append(
        Variable(
            "ncount",
            meta={
                "field_id": "ncount",
                "scale": "scaleLinear",
                "name": "N count",
                "datatype": "integer",
                "range": [min(values["n_counts"]), max(values["n_counts"])],
            },
            values=values["n_counts"],
            parents=[],
        )
    )
    if "z" not in kwargs["meta"].plot:
        kwargs["meta"].plot.update({"z": "length"})
    if "x" not in kwargs["meta"].plot:
        kwargs["meta"].plot.update({"x": "gc"})
    cov_range = [math.inf, -math.inf]
    read_cov_range = [math.inf, -math.inf]
    for cov_lib, cov_meta in blob_db["covLibs"].items():
        cov_file_name = field_name_from_path(blob_db["covLibs"][cov_lib]["f"])
        covs = values["%s_cov" % cov_lib]
        read_covs = values["%s_read_cov" % cov_lib]
        cov_range = [min(covs + [cov_range[0]]), max(covs + [cov_range[1]])]
        read_cov_range = [
            min(read_covs + [read_cov_range[0]]),
            max(read_covs + [read_cov_range[1]]),
        ]
        if "y" not in kwargs["meta"].plot:
            kwargs["meta"].plot.update({"y": "%s_cov" % cov_file_name})
        parsed.append(
            Variable(
                "%s_cov" % cov_file_name,
                values=covs,
                meta={"field_id": "%s_cov" % cov_file_name, "file": cov_meta["f"]},
                parents=cov.parent()
                + [
                    "children",
                    {
                        "id": "base_coverage",
                        "clamp": 1 if cov_range[0] == 0 else False,
                        "range": cov_range,
                    },
                    "children",
                ],
            )
        )
        parsed.append(
            Variable(
                "%s_read_cov" % cov_file_name,
                values=read_covs,
                meta={
                    "field_id": "%s_read_cov" % cov_file_name,
                    "file": cov_meta["f"],
                    "reads_mapped": cov_meta["reads_mapped"],
                    "reads_unmapped": cov_meta["reads_unmapped"],
                },
                parents=cov.parent()
                + [
                    "children",
                    {
                        "id": "read_coverage",
                        "datatype": "integer",
                        "clamp": 1 if read_cov_range[0] == 0 else False,
                        "range": read_cov_range,
                    },
                    "children",
                ],
            )
        )
    ranks = blob_db["dict_of_blobs"][identifiers.values[0]]["taxonomy"][
        blob_db["taxrules"][0]
    ].keys()
    for tax_rule in blob_db["taxrules"]:
        if "cat" not in kwargs["meta"].plot:
            kwargs["meta"].plot.update({"cat": "%s_phylum" % tax_rule})
        hit_list = hits_from_blob_db(blob_db, tax_rule)
        parsed.append(
            MultiArray(
                "%s_hits" % tax_rule,
                values=hit_list,
                meta={
                    "field_id": "%s_hits" % tax_rule,
                    "type": "multiarray",
                    "datatype": "mixed",
                    "preload": False,
                    "active": False,
                    "files": [m["f"] for x, m in blob_db["hitLibs"].items()],
                },
                parents=hits.parent() + ["children", {"id": tax_rule}, "children"],
                category_slot=None,
                headers=["taxid", "score"],
            )
        )
        for rank in ranks:
            field_id = "%s_%s" % (tax_rule, rank)
            parsed.append(
                Category(
                    field_id,
                    values=values[field_id],
                    meta={"field_id": field_id},
                    parents=hits.parent() + ["children", {"id": tax_rule}, "children"],
                )
            )
            parents = hits.parent() + [
                "children",
                {"id": tax_rule},
                "children",
                {"id": field_id},
                "data",
            ]
            field_id = "%s_%s_cindex" % (tax_rule, rank)
            parsed.append(
                Variable(
                    field_id,
                    values=values[field_id],
                    meta={
                        "scale": "scaleLinear",
                        "field_id": field_id,
                        "datatype": "integer",
                        "range": [min(values[field_id]), max(values[field_id])],
                        "preload": False,
                        "active": False,
                    },
                    parents=parents,
                )
            )
            field_id = "%s_%s_score" % (tax_rule, rank)
            _min = min(values[field_id])
            parsed.append(
                Variable(
                    field_id,
                    values=values[field_id],
                    meta={
                        "scale": "scaleLog",
                        "field_id": field_id,
                        "clamp": 1 if _min == 0 else False,
                        "datatype": "float",
                        "range": [_min, max(values[field_id])],
                        "preload": False,
                        "active": False,
                    },
                    parents=parents,
                )
            )

    return parsed


def parent():
    """Set standard metadata for synonyms."""
    return []
