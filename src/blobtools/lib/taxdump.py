#!/usr/bin/env python3
"""Taxdump Class module."""

import os
import re

from ..lib import file_io


class Taxdump:
    """Class for working with NCBI taxonomy."""

    __slots__ = ["directory", "ancestors", "names", "ranks"]

    def __init__(self, directory, **kwargs):
        """Init Taxdump class."""
        self.directory = directory
        self.ancestors = {}
        self.ranks = {}
        self.names = {}
        if kwargs:
            self.update_data(**kwargs)
        else:
            self.load_ranks()
            self.load_names()
            self.load_ancestors()

    def update_data(self, **kwargs):
        """Update values and keys for an existing field."""
        for key, value in kwargs.items():
            setattr(self, key, {int(taxid): data for taxid, data in value.items()})

    def load_ranks(self):
        """Load ranks from file."""
        filename = os.path.abspath(os.path.join(self.directory, "nodes.dmp"))
        try:
            for line in file_io.stream_file(filename):
                row = self.parse_taxdump_row(line)
                if len(row) > 1:
                    self.ranks[int(row[0])] = row[2]
        except TypeError:
            print("ERROR: Unable to parse %s." % filename)
            exit(1)

    def load_names(self):
        """Load names from file."""
        filename = os.path.abspath(os.path.join(self.directory, "names.dmp"))
        for line in file_io.stream_file(filename):
            row = self.parse_taxdump_row(line)
            if row and row[3] == "scientific name":
                self.names[int(row[0])] = row[1]

    def load_ancestors(self):
        """Load ancestors from file."""
        filename = os.path.abspath(os.path.join(self.directory, "taxidlineage.dmp"))
        for line in file_io.stream_file(filename):
            row = self.parse_taxdump_row(line)
            if row[1]:
                taxid = int(row[0])
                self.ancestors[taxid] = {
                    self.ranks[int(id)]: int(id)
                    for id in row[1].split(" ")
                    if self.ranks[int(id)] in self.list_ranks()
                }
                if self.ranks[taxid] in self.list_ranks():
                    self.ancestors[taxid].update({self.ranks[taxid]: taxid})
                last = 0
                for rank in self.list_ranks():
                    if rank in self.ancestors[taxid]:
                        last = -self.ancestors[taxid][rank]
                    else:
                        self.ancestors[taxid].update({rank: last})

    def values_to_dict(self):
        """Create a dict of values."""
        data = {}
        for key in ("ancestors", "names", "ranks"):
            if hasattr(self, key):
                data[key] = getattr(self, key)
        return data

    def lineage(self, taxid):
        """Create a dict of rank-name pairs for a taxid."""
        lineages = {}
        try:
            ancestors = self.ancestors[taxid]
        except ValueError:
            return {}
        for rank in self.list_ranks():
            if rank in ancestors and ancestors[rank] > 0:
                lineages.update({rank: self.names[ancestors[rank]]})
        return lineages

    @staticmethod
    def parse_taxdump_row(line):
        """Parse an ncbi taxdump file."""
        return re.split(r"\s*\|\s*", line)[:-1]

    @staticmethod
    def list_ranks():
        """Return a list of taxonomic ranks."""
        return [
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]


if __name__ == "__main__":
    import doctest

    doctest.testmod()
