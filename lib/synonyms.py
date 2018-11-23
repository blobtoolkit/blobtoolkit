#!/usr/bin/env python3
"""Parse BUSCO results into MultiArray Field."""

from pathlib import Path
import file_io
from field import Array


def parse(synonym_file, identifiers, **_kwargs):
    """Parse synonyms into Array."""
    meta = {}
    file_stem = Path(synonym_file).stem
    meta['field_id'] = "%s_synonyms" % file_stem
    by_id = {}
    ids = identifiers.to_set()
    for line in file_io.stream_file(synonym_file):
        row = line.rstrip().replace('"', '').split('\t')
        key = None
        names = []
        for value in row:
            if value in ids:
                key = value
            else:
                names.append(value)
        by_id.update({key: names})
    values = [by_id[id] if id in by_id else [] for id in identifiers.values]
    synonyms_field = Array(meta['field_id'], meta=meta, values=values, parents=['children'])
    return synonyms_field


def parent():
    """Set standard metadata for synonyms."""
    synonyms = {
        'datatype': 'string',
        'type': 'array',
        'id': 'synonyms',
        'name': 'Synonyms'
    }
    return [
        synonyms
    ]
