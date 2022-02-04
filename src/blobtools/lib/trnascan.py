#!/usr/bin/env python3
"""Parse BUSCO results into MultiArray Field."""

import re
from collections import defaultdict

from ..lib import file_io
from .field import MultiArray


def parse_trnascan(trnascan_file, identifiers):
    """Parse tRNAscan results into a MultiArray."""
    data = file_io.read_file(trnascan_file)
    lines = data.split("\n")
    header = True
    meta = {"file": trnascan_file}
    results = defaultdict(list)
    for line in lines:
        if header:
            row = re.split(" +", line)
            if len(row) > 1:
                if row[1].startswith("v."):
                    meta.update({"version": row[1]})
                elif row[1] == "Mode:":
                    meta.update({"mode": row[2]})
                    meta.update({"field_id": "trnascan_%s" % row[2].lower()})
                elif row[1].startswith("------"):
                    header = False
        else:
            row = re.split(r" +|\t", line)
            if len(row) == 9:
                results[row[0]].append([row[4], row[5]])
    if not identifiers.validate_list(list(results.keys())):
        raise UserWarning(
            "Contig names in the tRNAScan file did not match dataset identifiers."
        )
    values = [results[id] if id in results else [] for id in identifiers.values]
    trnascan_field = MultiArray(
        meta["field_id"],
        values=values,
        meta=meta,
        headers=("tRNA_type", "Anticodon"),
        parents=["children"],
    )
    return trnascan_field


def parse(files, **kwargs):
    """Parse all tRNAscan files."""
    parsed = []
    for file in files:
        parsed.append(
            parse_trnascan(file, identifiers=kwargs["dependencies"]["identifiers"])
        )
    return parsed


def parent():
    """Set standard metadata for tRNAscan."""
    trnascan = {
        "datatype": "mixed",
        "type": "array",
        "id": "trnascan",
        "name": "tRNAscan",
    }
    return [trnascan]
