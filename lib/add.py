#!/usr/bin/env python3

# pylint: disable=no-member

"""
Add data to a BlobDir.

Usage:
    blobtools add [--busco TSV...] [--hits TSV...] [--key path=value...]
                  [--link path=url...]
                  [--synonyms TSV...]
                  [--taxdump DIRECTORY] [--taxrule bestsum|bestsumorder]
                  [--replace] DIRECTORY

Arguments:
    DIRECTORY             Existing Blob directory.

Options:
    --busco TSV           BUSCO full_table.tsv output file.
    --hits TSV            Tabular BLAST/Diamond output file.
    --key path=value      Set a metadata key to value.
    --link path=url       Link to an external resource.
    --synonyms TSV        TSV file containing current identifiers and synonyms.
    --taxdump DIRECTORY   Location of NCBI new_taxdump directory.
    --taxrule bestsum|bestsumorder
                          Rule to use when assigning BLAST hits to taxa. [Default: bestsum]
    --replace             Replace existing fields with matching ids.

Examples:
    # 1. Add BUSCO scores to BlobDir
    ./blobtools add --busco busco.full_table.tsv BlobDir

"""

from docopt import docopt
import file_io
import hits
import busco
import key
import link
import synonyms
from taxdump import Taxdump
from field import Identifier
from dataset import Metadata

FIELDS = [{'flag': '--hits', 'module': hits},
          {'flag': '--busco', 'module': busco},
          {'flag': '--synonyms', 'module': synonyms}]
PARAMS = set(['--taxrule'])


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
    taxdump = None
    for field in FIELDS:
        for file in args[field['flag']]:
            if field['flag'] == '--hits':
                if not taxdump:
                    taxdump = fetch_taxdump(args['--taxdump'])
            parents = field['module'].parent()
            parsed = field['module'].parse(
                file,
                identifiers,
                **{key: args[key] for key in PARAMS},
                taxdump=taxdump)
            if not isinstance(parsed, list):
                parsed = [parsed]
            for data in parsed:
                if not args['--replace']:
                    if has_field_warning(meta, data.field_id):
                        continue
                meta.add_field(parents+data.parents, **data.meta)
                json_file = "%s/%s.json" % ('tests/files/dataset', data.field_id)
                file_io.write_file(json_file, data.values_to_dict())
    for string in args['--link']:
        link.add(string, meta, identifiers.values)
    for string in args['--key']:
        key.add(string, meta)
    file_io.write_file("%s/meta.update.json" % args['DIRECTORY'], meta.to_dict())


if __name__ == '__main__':
    main()
