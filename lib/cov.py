#!/usr/bin/env python3

# pylint: disable=no-member, too-many-locals

"""Parse base and read coverage into Variable Field."""

import os
import math
import time
from collections import defaultdict
from pathlib import Path
import pysam
from field import Variable


def parse_bam(bam_file, **kwargs):
    """Parse coverage into a Variables."""
    identifiers = kwargs['dependencies']['identifiers']
    lengths = kwargs['dependencies']['length'].values
    ncounts = kwargs['dependencies']['ncount'].values
    base_name = Path(bam_file).stem.split('.')[-1]
    index_file = Path("%s.bai" % bam_file)
    if not index_file.is_file():
        pysam.index(bam_file)
    else:
        index_file = False
    samfile = pysam.AlignmentFile(bam_file, "rb")
    # for read in samfile.fetch(identifiers.values[0], 100, 120):
    #     print(read)
    # samfile.close()
    _covs = defaultdict(int)
    _read_covs = {}
    ctr = 0
    start = time.time()
    for seq_id in identifiers.values:
        reads = set()
        for pileupcolumn in samfile.pileup(seq_id):
            _covs[seq_id] += pileupcolumn.n
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    reads.add(pileupread.alignment.query_name)
        _read_covs[seq_id] = len(reads)
        ctr += 1
        if ctr == 200:
            break
    print("%.1fs elapsed parsing %d contigs" % (time.time() - start, ctr))
    samfile.close()
    # stats = pysam.flagstat(bam_file)
    # print(stats)
    if index_file:
        os.remove(index_file)
    if not identifiers.validate_list(list(_covs.keys())):
        raise UserWarning('Contig names in the coverage file did not match dataset identifiers.')
    covs = []
    read_covs = []
    for index, seq_id in enumerate(identifiers.values):
        acgt_count = lengths[index] - ncounts[index]
        covs.append(float("%.4f" % (_covs[seq_id]/acgt_count)) if seq_id in _covs else 0)
        read_covs.append(_read_covs[seq_id] if seq_id in _read_covs else 0)
    field_id = "%s_cov" % base_name
    fields = {}
    fields['cov_range'] = [min(covs+[kwargs['cov_range'][0]]),
                           max(covs+[kwargs['cov_range'][1]])]
    fields['cov'] = Variable(field_id,
                             values=covs,
                             meta={'field_id': field_id},
                             parents=['children',
                                      {'id': 'base_coverage',
                                       'range': fields['cov_range']},
                                      'children']
                             )
    field_id = "%s_read_cov" % base_name
    fields['read_cov_range'] = [min(read_covs+[kwargs['read_cov_range'][0]]),
                                max(read_covs+[kwargs['read_cov_range'][1]])]
    fields['read_cov'] = Variable(field_id,
                                  values=read_covs,
                                  meta={'field_id': field_id},
                                  parents=['children',
                                           {'id': 'read_coverage',
                                            'type': 'integer',
                                            'range': fields['read_cov_range']},
                                           'children']
                                  )
    return fields


def parse(files, **kwargs):
    """Parse all BAM files."""
    parsed = []
    cov_range = [math.inf, -math.inf]
    read_cov_range = [math.inf, -math.inf]
    for file in files:
        fields = parse_bam(file, **kwargs, cov_range=cov_range, read_cov_range=read_cov_range)
        parsed.append(fields['cov'])
        parsed.append(fields['read_cov'])
        cov_range = fields['cov_range']
        read_cov_range = fields['read_cov_range']
    return parsed


def parent():
    """Set standard metadata for Coverage."""
    coverage = {
        'datatype': 'variable',
        'type': 'float',
        'scale': 'scaleLog',
        'id': 'coverage',
        'name': 'coverage'
    }
    return [
        coverage
    ]
