#!/usr/bin/env python3
"""Parse BUSCO results into MultiArray Field."""

import re
from collections import defaultdict
import file_io
from field import MultiArray


def parse_busco(busco_file, identifiers):
    """Parse BUSCO results into a MultiArray."""
    data = file_io.read_file(busco_file)
    lines = data.split('\n')
    rows = [re.split('\t', line) for line in lines[5:]]
    meta = {
        'version': lines[0].split(':')[1].strip(),
        'set': re.split(r':|\(|\)', lines[1])[1].strip(),
        'count': int(re.split(r':|\(|\)', lines[1])[5].strip()),
        'command': lines[2].split(':')[1].strip(),
        'file': busco_file
    }
    meta['set'] = re.search('-l.+\/(\w+_odb\d+)\/',meta['command'])[1]
    meta['field_id'] = "%s_busco" % meta['set']
    columns = re.split(r'# |\t', lines[4])[1:]
    busco_index = columns.index('Busco id')
    status_index = columns.index('Status')
    contig_index = columns.index('Contig')
    results = defaultdict(list)
    for row in rows:
        if len(row) > contig_index:
            results[row[contig_index]].append([row[busco_index], row[status_index]])
    if not identifiers.validate_list(list(results.keys())):
        raise UserWarning('Contig names in the Busco file did not match dataset identifiers.')
    values = [results[id] if id in results else [] for id in identifiers.values]
    busco_field = MultiArray(meta['field_id'],
                             values=values,
                             meta=meta,
                             headers=('Busco id', 'Status'),
                             parents=['children'],
                             category_slot=1
                             )
    return busco_field


def parse(files, **kwargs):
    """Parse all BUSCO files."""
    parsed = []
    for file in files:
        parsed.append(parse_busco(file, identifiers=kwargs['dependencies']['identifiers']))
    return parsed


def parent():
    """Set standard metadata for BUSCO."""
    busco = {
        'datatype': 'mixed',
        'type': 'array',
        'id': 'busco',
        'name': 'Busco'
    }
    return [
        busco
    ]
