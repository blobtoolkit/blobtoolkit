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
    meta['set'] = re.search(r'-l\s.*?\/*(\w+_odb\d+)\/', meta['command'])[1]
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


def busco_score(values, total):
    """Calculate BUSCO score."""
    fragmented = set()
    buscos = defaultdict(int)
    for contig in values:
        for busco in contig:
            if busco[1] == 'Fragmented':
                fragmented.add(busco[0])
            buscos[busco[0]] += 1
    present = len(buscos.keys())
    complete = present - len(fragmented)
    duplicated = len([value for value in buscos.values() if value > 1])
    scores = {
        'c': complete,
        'd': duplicated,
        'm': total - present,
        'f': len(fragmented),
        't': total,
        's': complete - duplicated
    }
    string = 'C:{:.1%}'.format(scores['c'] / total)
    string += '[S:{:.1%},'.format(scores['s'] / total)
    string += 'D:{:.1%}],'.format(scores['d'] / total)
    string += 'F:{:.1%},'.format(scores['f'] / total)
    string += 'M:{:.1%},'.format(scores['m'] / total)
    string += 'n:{:d}'.format(total)
    scores.update({'string': string})
    return scores


def summarise(indices, fields, **kwargs):  # pylint: disable=unused-argument
    """Summarise BUSCOs."""
    summary = {}
    for lineage in fields['lineages']:
        values = fields[lineage].expand_values()
        values = [values[i] for i in indices]
        total = fields[lineage].meta['count']
        lineage = lineage.replace('_busco', '')
        summary.update({lineage: busco_score(values, total)})
    return summary
