#!/usr/bin/env python3
"""Parse BLAST results into MultiArray Field."""

from collections import defaultdict
import file_io
from field import Category, MultiArray, Variable


def parse_blast(blast_file):
    """Parse file into dict of lists."""
    results = defaultdict(list)
    for line in file_io.stream_file(blast_file):
        row = line.rstrip().split('\t')
        seq_id, *offset = row[0].split('_-_')
        offset = int(offset[0]) if offset else 0
        hit = {'score': int(row[2]), 'start': int(row[9])+offset, 'end': int(row[10])+offset}
        try:
            hit.update({'taxid': int(row[1])})
        except ValueError:
            hit.update({'taxid': 0})
        results[seq_id].append(hit)
    return results


def apply_taxrule(blast, taxdump, taxrule, identifiers):
    """Apply taxrule to parsed BLAST results."""
    # results = [
    #     {'field_id': "%s_%s" % (taxrule, rank), 'values': defaultdict(str), 'data': {
    #         'cindex': defaultdict(int), 'score': defaultdict(int), 'positions': defaultdict(list)
    #     }}
    #     for rank in taxdump.list_ranks()
    # ]
    results = [
        {'field_id': "%s_%s" % (taxrule, rank), 'values': [], 'data': {
            'cindex': [], 'score': [], 'positions': []
        }}
        for rank in taxdump.list_ranks()
    ]
    values = [{
        'category': defaultdict(str),
        'cindex': defaultdict(int),
        'score': defaultdict(int),
        'positions': defaultdict(list)
        }
              for rank in taxdump.list_ranks()
             ]
    for seq_id, hits in blast.items():
        sorted_hits = sorted(hits, key=lambda k: k['score'], reverse=True)
        for index, rank in enumerate(taxdump.list_ranks()):
            cat_scores = defaultdict(int)
            for hit in sorted_hits:
                try:
                    category = taxdump.ancestors[hit['taxid']][rank]
                except KeyError:
                    category = 0
                if category > 0:
                    category = taxdump.names[category]
                elif category < 0:
                    category = "%s-undef" % taxdump.names[-category]
                else:
                    category = 'no-hit'
                values[index]['positions'][seq_id].append(
                    # results[index]['data']['positions'][seq_id].append(
                    [category, hit['start'], hit['end'], hit['score']]
                )
                if len(values[index]['positions'][seq_id]) < 10:
                    cat_scores[category] += hit['score']
            top_cat = max(cat_scores, key=cat_scores.get)
            values[index]['category'][seq_id] = top_cat
            values[index]['score'][seq_id] = cat_scores.get(top_cat)
            values[index]['cindex'][seq_id] = len(cat_scores.keys()) - 1
    for index, rank in enumerate(taxdump.list_ranks()):
        if not identifiers.validate_list(list(values[index]['category'].keys())):
            raise UserWarning('Contig names in the hits file do not match dataset identifiers.')
        results[index]['values'] = [
            values[index]['category'][id] if id in values[index]['category'] else 'no-hit'
            for id in identifiers.values]
        results[index]['data']['score'] = [
            values[index]['score'][id] if id in values[index]['score'] else 0
            for id in identifiers.values]
        results[index]['data']['cindex'] = [
            values[index]['cindex'][id] if id in values[index]['cindex'] else 0
            for id in identifiers.values]
        results[index]['data']['positions'] = [
            values[index]['positions'][id] if id in values[index]['positions'] else []
            for id in identifiers.values]
    return results


def parse(blast_file, identifiers, **kwargs):
    """Parse BLAST results into Fields."""
    blast = parse_blast(blast_file)
    results = apply_taxrule(blast, kwargs['taxdump'], kwargs['--taxrule'], identifiers)
    fields = []
    for result in results:
        parents = ['children']
        main = Category(result['field_id'],
                        values=result['values'],
                        meta={'field_id': result['field_id']},
                        parents=parents)
        fields.append(main)
        parents += [{'id': result['field_id']}, 'data']
        for subfield in ['score', 'cindex']:
            field_id = "%s_%s" % (result['field_id'], subfield)
            fields.append(Variable(field_id,
                                   values=result['data'][subfield],
                                   meta={'field_id': field_id},
                                   parents=parents))
        subfield = 'positions'
        field_id = "%s_%s" % (result['field_id'], subfield)
        fields.append(MultiArray(field_id,
                                 values=result['data'][subfield],
                                 fixed_keys=main.keys,
                                 meta={'field_id': field_id},
                                 parents=parents,
                                 category_slot=0,
                                 headers=['name', 'start', 'end', 'score']))
    return fields


def parent():
    """Set standard metadata for BLAST."""
    blast = {
        'datatype': 'mixed',
        'type': 'array',
        'id': 'blast',
        'name': 'blast',
        'children': []
    }
    return [
        blast
    ]
