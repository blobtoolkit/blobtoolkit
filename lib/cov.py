#!/usr/bin/env python3

# pylint: disable=no-member, too-many-locals

"""Parse base and read coverage into Variable Field."""

import os
import math
from multiprocessing import Pool
from pathlib import Path
import pysam
from tqdm import tqdm
from field import Variable
from run_external import seqtk_subseq


def _get_coverage(args):
    bam_file = args[0]
    f_char = args[1]
    seq_id = args[2]
    samfile = pysam.AlignmentFile(bam_file, "r%s" % f_char)
    cov = 0
    reads = set()
    for pileupcolumn in samfile.pileup(seq_id):
        cov += pileupcolumn.n
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                reads.add(pileupread.alignment.query_name)
    samfile.close()
    return [seq_id, cov, len(reads)]


def parse_bam(bam_file, **kwargs):
    """Parse coverage into a Variables."""
    identifiers = kwargs['dependencies']['identifiers']
    ids = identifiers.values
    lengths = kwargs['dependencies']['length'].values
    ncounts = kwargs['dependencies']['ncount'].values
    base_name = Path(bam_file).stem.split('.')[-1]
    f_char = Path(bam_file).suffix[1]
    index_file = Path("%s.bai" % bam_file)
    if not index_file.is_file():
        pysam.index(bam_file)
    else:
        index_file = False
    samfile = pysam.AlignmentFile(bam_file, "r%s" % f_char)
    stats = {'mapped': samfile.mapped,
             'unmapped': samfile.unmapped}
    samfile.close()
    print(stats)
    # samfile = pysam.AlignmentFile(bam_file, "r%s" % filetype_letter)
    print("Loading mapping data from %s" % bam_file)
    with Pool(int(kwargs['--threads'])) as pool:
        results = list(tqdm(pool.imap(_get_coverage,
                                      map(lambda x: (bam_file, f_char, x), ids)),
                            total=len(ids)))
    _covs = {}
    _read_covs = {}
    for result in results:
        _covs.update({result[0]: result[1]})
        _read_covs.update({result[0]: result[2]})
    if index_file:
        os.remove(index_file)
    if not identifiers.validate_list(list(_covs.keys())):
        raise UserWarning('Contig names in the coverage file did not match dataset identifiers.')
    covs = []
    read_covs = []
    for index, seq_id in enumerate(ids):
        acgt_count = lengths[index] - ncounts[index]
        covs.append(float("%.4f" % (_covs[seq_id]/acgt_count)) if seq_id in _covs else 0)
        read_covs.append(_read_covs[seq_id] if seq_id in _read_covs else 0)
    field_id = "%s_cov" % base_name
    fields = {}
    fields['cov_range'] = [min(covs+[kwargs['cov_range'][0]]),
                           max(covs+[kwargs['cov_range'][1]])]
    fields['cov'] = Variable(field_id,
                             values=covs,
                             meta={'field_id': field_id, 'file': bam_file},
                             parents=['children',
                                      {'id': 'base_coverage',
                                       'clamp': 1 if fields['cov_range'][0] == 0 else False,
                                       'range': fields['cov_range']},
                                      'children']
                             )
    field_id = "%s_read_cov" % base_name
    fields['read_cov_range'] = [min(read_covs+[kwargs['read_cov_range'][0]]),
                                max(read_covs+[kwargs['read_cov_range'][1]])]
    fields['read_cov'] = Variable(field_id,
                                  values=read_covs,
                                  meta={'field_id': field_id,
                                        'file': bam_file,
                                        'reads_total': stats['mapped'] + stats['unmapped'],
                                        'reads_mapped': stats['mapped'],
                                        'reads_unmapped': stats['unmapped']
                                        },
                                  parents=['children',
                                           {'id': 'read_coverage',
                                            'datatype': 'integer',
                                            'clamp': 1 if fields['read_cov_range'][0] == 0 else False,
                                            'range': fields['read_cov_range']},
                                           'children']
                                  )
    return fields


def apply_filter(ids, fastq_files, **kwargs):
    """Filter FASTQ file based on read alignment file."""
    suffix = kwargs['--suffix']
    bam_file = kwargs['--cov']
    read_ids = set()
    filetype_letter = Path(bam_file).suffix[1]
    index_file = Path("%s.bai" % bam_file)
    if not index_file.is_file():
        pysam.index(bam_file)
    else:
        index_file = False
    samfile = pysam.AlignmentFile(bam_file, "r%s" % filetype_letter)
    for seq_id in tqdm(ids):
        for read in samfile.fetch(seq_id):
            read_ids.add(read.query_name)
            break
    samfile.close()
    if index_file:
        os.remove(index_file)
    # read_ids = '\n'.join(list(read_ids))
    for fastq_file in fastq_files:
        path = Path(fastq_file)
        outfile = path.parent / (path.stem + '.' + suffix + path.suffix)
        seqtk_subseq(fastq_file, '\n'.join(list(read_ids)), outfile)


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
        'datatype': 'float',
        'type': 'variable',
        'scale': 'scaleLog',
        'id': 'coverage',
        'name': 'coverage'
    }
    return [
        coverage
    ]
