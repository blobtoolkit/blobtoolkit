#!/usr/bin/env python3

# pylint: disable=no-member

"""
Add data to a BlobDir.

Usage:
    blobtools add [--busco TSV...] [--hits TSV...]  [--fasta FASTA]
                  [--key path=value...] [--link path=url...] [--skip-link-test]
                  [--meta YAML] [--synonyms TSV...]
                  [--taxdump DIRECTORY] [--taxrule bestsum|bestsumorder]
                  [--create] [--replace] DIRECTORY

Arguments:
    DIRECTORY             Existing Blob directory.

Options:
    --busco TSV           BUSCO full_table.tsv output file.
    --fasta FASTA         FASTA sequence file.
    --hits TSV            Tabular BLAST/Diamond output file.
    --key path=value      Set a metadata key to value.
    --link path=URL       Link to an external resource.
    --skip-link-test      Skip test to see if link URL can be resolved.
    --meta YAML           Dataset metadata.
    --synonyms TSV        TSV file containing current identifiers and synonyms.
    --taxdump DIRECTORY   Location of NCBI new_taxdump directory.
    --taxrule bestsum|bestsumorder
                          Rule to use when assigning BLAST hits to taxa. [Default: bestsum]
    --create              Create a new BlobDir.
    --replace             Replace existing fields with matching ids.

Examples:
    # 1. Add BUSCO scores to BlobDir
    ./blobtools add --busco busco.full_table.tsv BlobDir

"""

import os
import glob
from docopt import docopt
import file_io
import busco
import fasta
import hits
import key
import link
import synonyms
from taxdump import Taxdump
from field import Identifier
from dataset import Metadata

FIELDS = [{'flag': '--fasta', 'module': fasta},
          {'flag': '--busco', 'module': busco},
          {'flag': '--hits', 'module': hits},
          {'flag': '--synonyms', 'module': synonyms}]
PARAMS = set(['--taxrule'])


def fetch_identifiers(path_to_dataset):
    """
    Load Identifiers from file.

    fetch_identifiers('tests/files/dataset')
    """
    try:
        data = file_io.load_yaml("%s/identifiers.json" % path_to_dataset)
        identifiers = Identifier('identifiers', **data)
    except TypeError:
        identifiers = False
    return identifiers


def fetch_metadata(path_to_dataset, **kwargs):
    """
    Load Metadata from file.

    fetch_metadata('tests/files/dataset')
    """
    dataset_id = path_to_dataset.split('/').pop()
    if kwargs['--meta']:
        meta = file_io.load_yaml(kwargs['--meta'])
        if kwargs['--replace']:
            files = glob.glob("%s/*" % kwargs['DIRECTORY'])
            for file in files:
                os.remove(file)
    else:
        meta = file_io.load_yaml("%s/meta.json" % path_to_dataset)
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


# def create_or_replace_dataset(directory, exists, create, replace):
#     """Warn about dataset existence."""
#     if exists:
#         if replace:
#             print("INFO: Replacing existing dataset \'%s\'." % directory)
#             return 'replace'
#         if create:
#             print("WARN: Dataset \'%s\' is already present, not overwriting." % directory)
#             print("WARN: Use '--replace' flag to overwrite existing dataset.")
#             exit()
#     if create:
#         print("INFO: Creating dataset \'%s\'." % directory)
#         return 'create'
#     print("WARN: Dataset \'%s\' does not exist." % directory)
#     print("WARN: Use '--create' flag to create a new dataset.")
#     exit()


def main():
    """Entrypoint for blobtools add."""
    args = docopt(__doc__)
    meta = fetch_metadata(args['DIRECTORY'], **args)
    identifiers = False
    taxdump = None
    for field in FIELDS:
        if args[field['flag']]:
            if not identifiers:
                identifiers = fetch_identifiers(args['DIRECTORY'])
            if field['flag'] == '--hits':
                if not taxdump:
                    taxdump = fetch_taxdump(args['--taxdump'])
            parents = field['module'].parent()
            parsed = field['module'].parse(
                args[field['flag']],
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
                json_file = "%s/%s.json" % (args['DIRECTORY'], data.field_id)
                file_io.write_file(json_file, data.values_to_dict())
    if not identifiers:
        identifiers = fetch_identifiers(args['DIRECTORY'])
    for string in args['--link']:
        link.add(string, meta, identifiers.values, args['--skip-link-test'])
    for string in args['--key']:
        key.add(string, meta)
    file_io.write_file("%s/meta.json" % args['DIRECTORY'], meta.to_dict())


if __name__ == '__main__':
    main()
