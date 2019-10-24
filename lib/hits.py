#!/usr/bin/env python3

# pylint: disable=too-many-locals,too-many-branches, unused-argument, unused-variable

"""Parse BLAST results into MultiArray Field."""

import math
from collections import defaultdict, Counter
import file_io
from field import Category, MultiArray, Variable


def parse_blast(blast_file, results=None, index=0):
    """Parse file into dict of lists."""
    if results is None:
        results = defaultdict(list)
    for line in file_io.stream_file(blast_file):
        row = line.rstrip().split('\t')
        seq_id, *offset = row[0].split('_-_')
        offset = int(offset[0]) if offset else 0
        try:
            hit = {'subject': row[4],
                   'score': float(row[2]),
                   'start': int(row[9])+offset,
                   'end': int(row[10])+offset,
                   'file': index}
        except IndexError:
            hit = {'subject': row[3],
                   'score': float(row[2]),
                   'start': None,
                   'end': None,
                   'file': index}
        try:
            hit.update({'taxid': int(row[1])})
        except ValueError:
            hit.update({'taxid': 0})
        results[seq_id].append(hit)
    return results


def apply_taxrule(blast, taxdump, taxrule, identifiers, results=None):
    """Apply taxrule to parsed BLAST results."""
    if results is None:
        blank = [None] * len(identifiers.values)
        results = [
            {'field_id': "%s_%s" % (taxrule, rank), 'values': blank[:], 'data': {
                'cindex': blank[:], 'score': blank[:], 'positions': blank[:],
                'hits': blank[:]
            }}
            for rank in taxdump.list_ranks()
        ]
    values = [{
        'category': defaultdict(str),
        'cindex': defaultdict(int),
        'score': defaultdict(float),
        'positions': defaultdict(list),
        'hits': defaultdict(list)
        } for rank in taxdump.list_ranks()]
    for seq_id, hits in blast.items():
        sorted_hits = sorted(hits, key=lambda k: k['score'], reverse=True)
        for index, rank in enumerate(taxdump.list_ranks()):
            cat_scores = defaultdict(float)
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
                    category = 'undef'
                values[index]['positions'][seq_id].append(
                    [category]
                )
                if index == 0:
                    try:
                        values[index]['hits'][seq_id].append(
                            [hit['taxid'], hit['start'], hit['end'], hit['score'], hit['subject'], hit['file']]
                            )
                    except KeyError:
                        values[index]['hits'][seq_id].append(
                            [category, hit['file']]
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
        for i, seq_id in enumerate(identifiers.values):
            if results[index]['data']['score'][i] in (0, None):
                if seq_id in values[index]['category']:
                    results[index]['values'][i] = values[index]['category'][seq_id]
                    results[index]['data']['score'][i] = values[index]['score'][seq_id]
                    results[index]['data']['cindex'][i] = values[index]['cindex'][seq_id]
                    results[index]['data']['positions'][i] = values[index]['positions'][seq_id]
                    if index == 0:
                        results[index]['data']['hits'][i] = values[index]['hits'][seq_id]
                else:
                    results[index]['values'][i] = 'no-hit'
                    results[index]['data']['score'][i] = 0
                    results[index]['data']['cindex'][i] = 0
                    results[index]['data']['positions'][i] = []
                    if index == 0:
                        results[index]['data']['hits'][i] = []
    return results


def create_fields(results, taxrule, files, fields=None):
    """Store BLAST results as Fields."""
    if fields is None:
        fields = []
    hits_id = "%s_%s" % (taxrule, 'positions')
    fields.append(MultiArray(hits_id,
                             values=results[0]['data']['hits'],
                             meta={
                                 'field_id': hits_id,
                                 'name': hits_id,
                                 'type': 'multiarray',
                                 'datatype': 'mixed',
                                 'preload': False,
                                 'active': False,
                                 'files': files
                                 },
                             parents=['children', {'id': taxrule}, 'children'],
                             category_slot=None,
                             headers=['taxid', 'start', 'end', 'score', 'subject', 'index']))
    for result in results:
        main = Category(result['field_id'],
                        values=result['values'],
                        meta={
                            'field_id': result['field_id'],
                            'name': result['field_id']
                        },
                        parents=['children', {'id': taxrule}, 'children'])
        fields.append(main)
        parents = ['children', {'id': taxrule}, 'children', {'id': result['field_id']}, 'data']
        field_id = "%s_%s" % (result['field_id'], 'cindex')
        fields.append(Variable(field_id,
                               values=result['data']['cindex'],
                               meta={
                                   'scale': 'scaleLinear',
                                   'field_id': field_id,
                                   'name': field_id,
                                   'datatype': 'integer',
                                   'range': [min(result['data']['cindex']),
                                             max(result['data']['cindex'])],
                                   'preload': False,
                                   'active': False
                                   },
                               parents=parents))
        field_id = "%s_%s" % (result['field_id'], 'score')
        _min = min(result['data']['score'])
        fields.append(Variable(field_id,
                               values=result['data']['score'],
                               meta={
                                   'scale': 'scaleLog',
                                   'field_id': field_id,
                                   'name': field_id,
                                   'clamp': 1 if _min == 0 else False,
                                   'datatype': 'float',
                                   'range': [_min,
                                             max(result['data']['score'])],
                                   'preload': False,
                                   'active': False
                                   },
                               parents=parents))
        subfield = 'positions'
        field_id = "%s_%s" % (result['field_id'], subfield)
        if len(result['data'][subfield]) > 1:
            headers = ['name']
        else:
            headers = ['name']
        fields.append(MultiArray(field_id,
                                 values=result['data'][subfield],
                                 fixed_keys=main.keys,
                                 meta={
                                     'field_id': field_id,
                                     'name': field_id,
                                     'type': 'multiarray',
                                     'datatype': 'string',
                                     'preload': False,
                                     'active': False,
                                     'linked_field': hits_id
                                     },
                                 parents=parents,
                                 category_slot=0,
                                 headers=headers))
    return fields


def parse(files, **kwargs):
    """Parse BLAST results into Fields."""
    blast = None
    fields = []
    identifiers = kwargs['dependencies']['identifiers']
    if kwargs['--taxrule'] == 'bestsum':
        for index, file in enumerate(files):
            blast = parse_blast(file, blast, index)
        results = apply_taxrule(blast, kwargs['taxdump'], kwargs['--taxrule'], identifiers)
        fields = create_fields(results, kwargs['--taxrule'], files)
    elif kwargs['--taxrule'] == 'bestsumorder':
        results = None
        for index, file in enumerate(files):
            blast = parse_blast(file, None, index)
            results = apply_taxrule(blast,
                                    kwargs['taxdump'],
                                    kwargs['--taxrule'],
                                    identifiers,
                                    results)
        fields = create_fields(results, kwargs['--taxrule'], files)
    if 'cat' not in kwargs['meta'].plot:
        kwargs['meta'].plot.update({'cat': "%s_phylum" % kwargs['--taxrule']})
    return fields


def parent():
    """Set standard metadata for BLAST."""
    blast = {
        'datatype': 'string',
        'type': 'category',
        'id': 'taxonomy',
        'name': 'Taxonomy'
    }
    return [
        blast
    ]


def summarise(indices, fields, **kwargs):
    """Summarise assembly sequence stats."""
    taxonomy = kwargs['stats']['taxonomy']
    summary = {}
    lengths = [fields['length'].values[i] for i in indices]
    gcs = [fields['gc'].values[i] for i in indices]
    covs = [fields['cov'].values[i] for i in indices] if 'cov' in fields else []
    summary.update({'total': length_stats(lengths, gcs, covs)})
    hits = list(map(lambda x: fields['hits'].keys[x], [fields['hits'].values[i] for i in indices]))
    counts = Counter(hits)
    index = 0
    limit = 9
    other = []
    other_gcs = []
    other_covs = []
    for key in sorted(counts.items(), key=lambda x: x[1], reverse=True):
        taxon = key[0]
        subset = [i for i, j in enumerate(hits) if j == taxon]
        lengths = [fields['length'].values[indices[i]] for i in subset]
        gcs = [fields['gc'].values[indices[i]] for i in subset]
        covs = [fields['cov'].values[indices[i]] for i in subset] if 'cov' in fields else []
        if len(counts.keys()) > limit and index >= limit:
            other += lengths
            other_gcs += gcs
            other_covs += covs
            if 'target' in taxonomy and taxon == taxonomy['target']:
                summary.update({'target': length_stats(lengths, gcs, covs)})
        else:
            summary.update({taxon: length_stats(lengths, gcs, covs)})
        index += 1
    if other:
        summary.update({'other': length_stats(other, other_gcs, other_covs)})
    return summary


def weighted_mean_sd(values, weights, log=False):
    """Calculate weighted mean, median and standard deviation."""
    weighted = []
    total_weight = 0
    total = 0
    for weight, value in zip(weights, values):
        if log:
            value = math.log10(value) + 10
        weighted.append({'weight': weight, 'value': value})
        total_weight += weight
        total += value * weight
    mid = total_weight / 2
    running_total = 0
    median = 0
    for entry in sorted(weighted, key=lambda i: i['value']):
        running_total += entry['weight']
        if running_total > mid:
            median = entry['value']
            break
    mean = total / total_weight
    if median == 0:
        median = mean
    variance = sum([entry['weight'] * ((entry['value'] - mean) ** 2)
                    for entry in weighted]
                   ) / total_weight
    st_dev = variance ** 0.5
    upper = mean + 2 * st_dev
    lower = mean - 2 * st_dev
    if log:
        upper = 10 ** (upper - 10)
        lower = 10 ** (lower - 10)
        mean = 10 ** (mean - 10)
        median = 10 ** (median - 10)
    return mean, median, st_dev, upper, lower


def length_stats(all_lengths, all_gcs, all_covs):
    """Calculate stats for a set of sequence lengths."""
    span = sum(all_lengths)
    count = len(all_lengths)
    lengths = []
    gcs = []
    covs = []
    if all_covs:
        for length, gc_value, cov in zip(all_lengths, all_gcs, all_covs):
            if cov >= 0.01:
                lengths.append(length)
                gcs.append(gc_value)
                covs.append(cov)
    else:
        lengths = all_lengths
        gcs = all_gcs
    gc_mean, gc_median, gc_dev, gc_upper, gc_lower = weighted_mean_sd(gcs, lengths)
    stats = {'span': span,
             'count': count,
             'gc': [float("%.4f" % gc_mean),
                    float("%.4f" % gc_median),
                    float("%.4f" % gc_lower),
                    float("%.4f" % gc_upper),
                    float("%.4f" % min(gcs)),
                    float("%.4f" % max(gcs))]
             }
    if covs:
        cov_mean, cov_median, cov_dev, cov_upper, cov_lower = weighted_mean_sd(covs, lengths, log=True)
        stats.update({'cov': [float("%.4f" % cov_mean),
                              float("%.4f" % cov_median),
                              float("%.4f" % cov_lower),
                              float("%.4f" % cov_upper),
                              float("%.4f" % min(covs)),
                              float("%.4f" % max(covs))]})
    n50 = span * 0.5
    n90 = span * 0.9
    all_lengths.sort(reverse=True)
    nlength = 0
    ncount = 0
    for length in all_lengths:
        ncount += 1
        nlength += length
        if 'n50' not in stats and nlength > n50:
            stats.update({'n50': length, 'l50': ncount})
        if 'n90' not in stats and nlength > n90:
            stats.update({'n90': length, 'l90': ncount})
    if 'n50' not in stats:
        stats.update({'n50': all_lengths[-1], 'l50': ncount})
    if 'n90' not in stats:
        stats.update({'n90': all_lengths[-1], 'l90': ncount})
    return stats
