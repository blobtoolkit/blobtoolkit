#!/usr/bin/env python3
"""Parse BUSCO results into MultiArray Field."""

import re
from pathlib import Path
import file_io
from field import Array


def parse(synonym_file, identifiers):
    """Parse synonyms into Array."""
    data = file_io.read_file(synonym_file)
    lines = data.split('\n')
    rows = [re.split('\t', line.replace('"', '')) for line in lines]
    meta = {}
    file_stem = Path(synonym_file).stem
    meta['field_id'] = "%s_synonyms" % file_stem
    by_id = {}
    ids = identifiers.to_set()
    for row in rows:
        key = None
        names = []
        for value in row:
            if value in ids:
                key = value
            else:
                names.append(value)
        by_id.update({key: names})
    values = [by_id[id] if id in by_id else [] for id in identifiers.values]
    synonyms_field = Array(meta['field_id'], meta=meta, values=values)
    return synonyms_field


def parent():
    """Set standard metadata for synonyms."""
    synonyms = {
        'datatype': 'string',
        'type': 'array',
        'id': 'synonyms',
        'name': 'Synonyms',
        'children': []
    }
    return [
        synonyms,
        'children'
    ]
