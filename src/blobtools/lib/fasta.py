#!/usr/bin/env python3
# pylint: skip-file

"""Parse FASTA sequence into Fields."""

from collections import Counter
from collections import OrderedDict

from blobtk import filter
from tqdm import tqdm

from ..lib import file_io
from .field import Identifier
from .field import Variable


def base_composition(seq_str):
    """Sequence base composition summary."""
    counts = Counter(seq_str.upper())
    at_bases = ["A", "T", "W"]
    gc_bases = ["C", "G", "S"]
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


def apply_filter(ids, fasta_file, **kwargs):
    """Filter FASTA format assembly."""
    suffix = kwargs["--suffix"]
    filter.fastx(list=ids, fasta=fasta_file, suffix=suffix, fasta_out=True)


def parse(file, **kwargs):
    """Parse all synonym files."""
    parsed = []
    _lengths = OrderedDict()
    _gc_portions = OrderedDict()
    _n_counts = OrderedDict()
    lengths = []
    gc_portions = []
    n_counts = []
    print(f"Loading sequences from {file}")
    pbar = tqdm(file_io.stream_fasta(file))
    for seq_id, seq_str in pbar:
        pbar.set_description(f" - processing {seq_id}")
        _lengths[seq_id] = len(seq_str)
        _gc_portions[seq_id], _n_counts[seq_id] = base_composition(seq_str)
    identifiers = kwargs["dependencies"]["identifiers"]
    if not identifiers:
        identifiers = Identifier(
            "identifiers",
            meta={"field_id": "identifiers"},
            values=list(_lengths.keys()),
            parents=[],
        )
        kwargs["meta"].assembly.update({"scaffold-count": len(identifiers.values)})
        parsed.append(identifiers)
    for seq_id in identifiers.values:
        lengths.append(_lengths[seq_id] if seq_id in _lengths else 0)
        gc_portions.append(_gc_portions[seq_id] if seq_id in _gc_portions else 0)
        n_counts.append(_n_counts[seq_id] if seq_id in _n_counts else 0)
    kwargs["meta"].assembly.update({"span": sum(lengths)})
    parsed.append(
        Variable(
            "gc",
            meta={
                "field_id": "gc",
                "preload": True,
                "scale": "scaleLinear",
                "name": "GC",
                "datatype": "float",
                "range": [min(gc_portions), max(gc_portions)],
            },
            values=gc_portions,
            parents=[],
        )
    )
    _min = min(lengths)
    parsed.append(
        Variable(
            "length",
            meta={
                "preload": True,
                "scale": "scaleLog",
                "field_id": "length",
                "name": "Length",
                "clamp": 1 if _min == 0 else False,
                "datatype": "integer",
                "range": [_min, max(lengths)],
            },
            values=lengths,
            parents=[],
        )
    )
    parsed.append(
        Variable(
            "ncount",
            meta={
                "scale": "scaleLinear",
                "field_id": "ncount",
                "name": "N count",
                "datatype": "integer",
                "range": [min(n_counts), max(n_counts)],
            },
            values=n_counts,
            parents=[],
        )
    )
    if "x" not in kwargs["meta"].plot:
        kwargs["meta"].plot.update({"x": "gc"})
    if "z" not in kwargs["meta"].plot:
        kwargs["meta"].plot.update({"z": "length"})
    return parsed


def parent():
    """Set standard metadata for synonyms."""
    return []


def summarise(indices, fields, **kwargs):
    """Summarise assembly sequence stats."""
    gcs = [fields["gc"].values[i] for i in indices]
    lengths = [fields["length"].values[i] for i in indices]
    gc_mean = sum([gc * length for gc, length in zip(gcs, lengths)]) / sum(lengths)
    gc = float("%.4f" % gc_mean)
    at = 1 - gc
    if "ncount" in fields and fields["ncount"]:
        ncount = sum([fields["ncount"].values[i] for i in indices])
    else:
        ncount = 0
    span = sum(lengths)
    return {"at": at, "gc": gc, "n": float("%.4f" % (ncount / span))}


def remove_from_meta(meta):
    """Delete all fasta-derived fields."""
    field_ids = []
    if meta.has_field("gc"):
        field_ids += meta.remove_field("gc")
        meta.plot.pop("x", None)
    if meta.has_field("length"):
        field_ids += meta.remove_field("length")
        meta.plot.pop("z", None)
    if meta.has_field("ncount"):
        field_ids += meta.remove_field("ncount")
    meta.assembly.pop("file", None)
    return field_ids
