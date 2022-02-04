#!/usr/bin/env python3

# pylint: disable=unused-argument, broad-except

"""Add taxonomic ranks to metadata for a taxid."""


def add(taxid, taxdump, meta):
    """Add ranks to meta."""
    taxid = int(taxid)
    lineage = taxdump.lineage(taxid)
    if lineage:
        meta.taxon.update({'taxid': taxid})
        for rank, value in lineage.items():
            meta.taxon.update({rank: value})
    else:
        print("WARN: '%s' was not found in the taxdump." % taxid)
    return meta.taxon


def summarise(indices, fields, **kwargs):
    """Summarise taxonomy."""
    summary = {}
    meta = kwargs['meta']
    try:
        if 'taxid' in meta.taxon:
            summary.update({'taxid': meta.taxon['taxid']})
        names = []
        for rank in ('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'):
            if rank in meta.taxon:
                names.append(meta.taxon[rank])
        if names:
            summary.update({'lineage': '; '.join(names)})
        if 'cat' in meta.plot:
            rank = meta.plot['cat'].split('_')[-1]
            summary.update({'targetRank': rank})
            summary.update({'target': meta.taxon[rank]})
    except Exception:
        pass
    return summary
