#!/usr/bin/env python3
"""Functions to fetch data and metadata."""

import os
import glob
import file_io
from taxdump import Taxdump
from field import Field, Identifier, Category, Variable, Array, MultiArray
from dataset import Metadata

TYPES = {
    'identifier': Identifier,
    'category': Category,
    'variable': Variable,
    'array': Array,
    'multiarray': MultiArray,
    'field': Field
}


def fetch_field(path_to_dataset, field_id):
    """
    Load fields from file.

    fetch_field('tests/files/dataset', 'identifiers')
    """
    try:
        data = file_io.load_yaml("%s/%s.json" % (path_to_dataset, field_id))
        field = TYPES[data['type']](field_id, **data)
    except TypeError:
        field = False
    except KeyError:
        field = False
    return field


def fetch_metadata(path_to_dataset, **kwargs):
    """
    Load Metadata from file.

    fetch_metadata('tests/files/dataset')
    """
    dataset_id = path_to_dataset.split('/').pop()
    if not os.path.exists(path_to_dataset):
        os.makedirs(path_to_dataset)
    if kwargs.get('--meta'):
        meta = file_io.load_yaml(kwargs['--meta'])
        if kwargs['--replace']:
            files = glob.glob("%s/*" % kwargs['DIRECTORY'])
            for file in files:
                os.remove(file)
    elif not kwargs.get('meta'):
        meta = file_io.load_yaml("%s/meta.json" % path_to_dataset)
    else:
        meta = kwargs['meta']
    if 'id' not in meta:
        meta['id'] = dataset_id
    return Metadata(dataset_id, **meta)


def fetch_taxdump(path_to_taxdump):
    """Load Taxdump from file."""
    json_file = "%s/taxdump.json" % path_to_taxdump
    data = file_io.load_yaml(json_file)
    if data is None:
        print('Parsing taxdump')
        taxdump = Taxdump(path_to_taxdump)
        file_io.write_file(json_file, taxdump.values_to_dict())
    else:
        print('Loading parsed taxdump')
        taxdump = Taxdump(path_to_taxdump, **data)
    return taxdump


if __name__ == '__main__':
    import doctest
    doctest.testmod()
