#!/usr/bin/env python3

"""
Add data to a BlobDir.

Usage:
    blobtools add [--busco TSV...] [--key path=value...] [--link path=url...]
                  [--synonyms TSV...] [--replace] DIRECTORY

Arguments:
    DIRECTORY             Existing Blob directory.

Options:
    --busco TSV           BUSCO full_table.tsv output file.
    --key path=value      Set a metadata key to value.
    --link path=url       Link to an external resource.
    --synonyms TSV        TSV file containing current identifiers and synonyms.
    --replace             Replace existing fields with matching ids.

Examples:
    # 1. Add BUSCO scores to BlobDir
    ./blobtools add --busco busco.full_table.tsv BlobDir

"""

from docopt import docopt
import file_io
import busco
import key
import link
import synonyms
from field import Identifier
from dataset import Metadata

FIELDS = [{'flag': '--busco', 'module': busco},
          {'flag': '--synonyms', 'module': synonyms}]


def fetch_identifiers(path_to_dataset):
    """
    Load Identifiers from file.

    fetch_identifiers('tests/files/dataset')
    """
    data = file_io.load_yaml("%s/identifiers.json" % path_to_dataset)
    return Identifier('identifiers', **data)


def fetch_metadata(path_to_dataset):
    """
    Load Metadata from file.

    fetch_metadata('tests/files/dataset')
    """
    meta = file_io.load_yaml("%s/meta.json" % path_to_dataset)
    dataset_id = path_to_dataset.split('/').pop()
    return Metadata(dataset_id, **meta)


def has_field_warning(meta, field_id):
    """Warn if dataset has existing field with same id."""
    if meta.has_field(field_id):
        print("WARN: Field \'%s\' is already present in dataset, not overwriting." % field_id)
        print("WARN: Use '--replace' flag to overwrite existing field.")
        return 1
    return 0


def main():
    """Entrypoint for blobtools add."""
    args = docopt(__doc__)
    meta = fetch_metadata(args['DIRECTORY'])
    identifiers = fetch_identifiers(args['DIRECTORY'])
    for field in FIELDS:
        for file in args[field['flag']]:
            data = field['module'].parse(file, identifiers)
            if not args['--replace']:
                if has_field_warning(meta, data.field_id):
                    continue
            parents = field['module'].parent()
            meta.add_field(parents, **data.meta)
            json_file = "%s/%s.json" % ('tests/files/dataset', data.field_id)
            file_io.write_file(json_file, data.values_to_dict())
    for string in args['--link']:
        link.add(string, meta, identifiers.values)
    for string in args['--key']:
        key.add(string, meta)
    file_io.write_file("%s/meta.update.json" % args['DIRECTORY'], meta.to_dict())


if __name__ == '__main__':
    main()
