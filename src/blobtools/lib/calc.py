#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-nested-blocks

"""
Field-based calculations.

Usage:
    blobtools calc [--combine STRING...] [--equation STRING...] [--sum LIST...]
                   [--replace] DIRECTORY

Arguments:
    DIRECTORY             Existing Blob directory.

Options:
    --combine   STRING    Group categories FIELD:CATEGORY,LIST:ALIAS.
    --equation  STRING    Combine variable fields using ()/*+- operators, e.g. base_cov*(length-ncount).
    --sum       LIST      Sum all fields (and/or child fields) in a comman separated list of field IDs.
    --replace             Replace existing fields with matching ids.

Examples:
    blobtools calc --sum base_coverage=summed_cov
    blobtools calc --equation ncount/length=n_fraction
    blobtools calc --equation ncount/length*100=n_percent
    blobtools calc --equation (base_cov1+base_cov2)*(length-ncount)=combined_base_count
    blobtools calc --combine bestdistorder_species:(bestdistorder_phylum:Nematoda):Nematoda=nematoda_grouped
    blobtools calc --combine bestdistorder_species:(bestdistorder_phylum:Nematoda,Arthropoda):Contaminants=contaminants
    blobtools calc --combine bestdistorder_genus:Homo,Pan:Apes=apes_grouped
    blobtools calc --combine bestdistorder_genus:Homo,Pan:Apes::bestdistorder_genus:Mus,Rattus:Rodents=apes_and_rodents
    blobtools calc --combine bestdistorder_genus:(bestdistorder_genus_cindex:>0):Uncertain=uncertain_grouped

"""

from docopt import docopt

# from fetch import fetch_metadata


def has_field_warning(meta, field_id):
    """Warn if dataset has existing field with same id."""
    if meta.has_field(field_id):
        print(
            "WARN: Field '%s' is already present in dataset, not overwriting."
            % field_id
        )
        print("WARN: Use '--replace' flag to overwrite existing field.")
        return 1
    return 0


def main(args):
    """Entrypoint for blobtools calc."""
    # meta = fetch_metadata(args['DIRECTORY'], **args)
    if args["--sum"]:
        print(args["--sum"])
    if args["--equation"]:
        print(args["--equation"])
    if args["--combine"]:
        print(args["--combine"])

    # if not args['--replace']:
    #     if has_field_warning(meta, data.field_id):
    #         continue
    #     for parent in data.parents:
    #         if 'range' in parent:
    #             parent_meta = meta.field_meta(parent['id'])
    #             if parent_meta and 'range' in parent_meta:
    #                 parent['range'][0] = min(parent['range'][0], parent_meta['range'][0])
    #                 parent['range'][1] = max(parent['range'][1], parent_meta['range'][1])
    #     meta.add_field(parents+data.parents, **data.meta)
    #     json_file = "%s/%s.json" % (args['DIRECTORY'], data.field_id)
    #     file_io.write_file(json_file, data.values_to_dict())
    # file_io.write_file("%s/meta.json" % args['DIRECTORY'], meta.to_dict())


if __name__ == "__main__":
    main()
