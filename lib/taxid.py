#!/usr/bin/env python3
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
    return True
