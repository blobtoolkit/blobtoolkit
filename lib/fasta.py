#!/usr/bin/env python3
"""Parse FASTA sequence into Fields."""

from collections import Counter, OrderedDict
import file_io
from field import Identifier, Variable


def base_composition(seq_str):
    """Sequence base composition summary."""
    counts = Counter(seq_str.upper())
    at_bases = ['A', 'T', 'W']
    gc_bases = ['C', 'G', 'S']
    at_count = 0
    gc_count = 0
    for base in at_bases:
        at_count += counts[base]
    for base in gc_bases:
        gc_count += counts[base]
    acgt_count = at_count + gc_count
    gc_portion = float("%.4f" % (gc_count / acgt_count))
    n_count = len(seq_str) - acgt_count
    return gc_portion, n_count


def parse(file, identifiers, **_kwargs):
    """Parse all synonym files."""
    parsed = []
    _lengths = OrderedDict()
    _gc_portions = OrderedDict()
    _n_counts = OrderedDict()
    lengths = []
    gc_portions = []
    n_counts = []
    for seq_id, seq_str in file_io.stream_fasta(file):
        _lengths[seq_id] = len(seq_str)
        _gc_portions[seq_id], _n_counts[seq_id] = base_composition(seq_str)
    if not identifiers:
        identifiers = Identifier('identifiers',
                                 meta={'field_id': 'identifiers'},
                                 values=list(_lengths.keys()),
                                 parents=[])
        parsed.append(identifiers)
    for seq_id in identifiers.values:
        lengths.append(_lengths[seq_id] if seq_id in _lengths else 0)
        gc_portions.append(_gc_portions[seq_id] if seq_id in _gc_portions else 0)
        n_counts.append(_n_counts[seq_id] if seq_id in _n_counts else 0)
    # else:
    #     lengths = list(_lengths.values())
    #     gc_portions = list(_gc_portions.values())
    #     n_counts = list(_n_counts.values())
    #     parsed.append(Identifier('identifiers',
    #                              meta={'field_id': 'identifiers'},
    #                              values=list(_lengths.keys()),
    #                              parents=[]))
    parsed.append(Variable('gc',
                           meta={
                               'preload': True,
                               'scale': 'scaleLinear',
                               'field_id': 'gc',
                               'name': 'GC',
                               'datatype': 'float',
                               'range': [min(gc_portions), max(gc_portions)]
                           },
                           values=gc_portions,
                           parents=[]))
    parsed.append(Variable('length',
                           meta={
                               'preload': True,
                               'scale': 'scaleLog',
                               'field_id': 'length',
                               'name': 'Length',
                               'datatype': 'integer',
                               'range': [min(lengths), max(lengths)]
                           },
                           values=lengths,
                           parents=[]))
    parsed.append(Variable('ncount',
                           meta={
                               'scale': 'scaleLinear',
                               'field_id': 'ncount',
                               'name': 'N count',
                               'datatype': 'integer',
                               'range': [min(n_counts), max(n_counts)]
                           },
                           values=n_counts,
                           parents=[]))

    return parsed


def parent():
    """Set standard metadata for synonyms."""
    return []
