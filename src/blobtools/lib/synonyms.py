#!/usr/bin/env python3
"""Parse BUSCO results into MultiArray Field."""

import re
from pathlib import Path

from ..lib import file_io
from .field import Array
from .text import parse_header_row
from .text import set_delimiter


def parse_synonyms(synonym_file, delimiter, columns, header, identifiers):
    """Parse synonyms into Array."""
    meta = {}
    synonym_file, *prefix = synonym_file.split("=")
    if prefix:
        prefix = prefix[0]
    else:
        prefix = Path(synonym_file).stem
    meta["field_id"] = "%s_synonyms" % prefix
    by_id = {}
    ids = identifiers.to_set()
    data = file_io.stream_file(synonym_file)
    lines = [line for line in data]
    if columns:
        columns = columns.split(",")
    else:
        columns = []
    delimit = set_delimiter(delimiter, sample=lines[0])
    if header:
        header_row = lines[0].rstrip().replace('"', "")
        columns = parse_header_row(delimit, header_row, columns)
        lines = lines[1:]
    try:
        id_col = columns.index("identifier")
    except ValueError:
        id_col = None
    for line in lines:
        row = re.split(delimit, line.rstrip().replace('"', ""))
        key = None
        names = []
        for i, value in enumerate(row):
            if id_col is None and value in ids:
                key = value
                id_col = i
            elif i == id_col:
                key = value
            else:
                names.append(value)
        by_id.update({key: names})
    values = [by_id[id] if id in by_id else [] for id in identifiers.values]
    del columns[id_col]
    synonyms_field = Array(
        meta["field_id"],
        meta=meta,
        values=values,
        headers=columns,
        parents=["children"],
    )
    return synonyms_field


def parse(files, **kwargs):
    """Parse all synonym files."""
    parsed = []
    for file in files:
        parsed.append(
            parse_synonyms(
                file,
                delimiter=kwargs["--text-delimiter"],
                columns=kwargs["--text-cols"],
                header=kwargs["--text-header"],
                identifiers=kwargs["dependencies"]["identifiers"],
            )
        )
    return parsed


def parent():
    """Set standard metadata for synonyms."""
    synonyms = {
        "datatype": "string",
        "type": "array",
        "id": "synonyms",
        "name": "Synonyms",
    }
    return [synonyms]
