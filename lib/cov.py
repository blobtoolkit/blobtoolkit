#!/usr/bin/env python3

# pylint: disable=no-member, too-many-locals
# pylint: skip-file

"""Parse base and read coverage into Variable Field."""

import os
import math
from multiprocessing import Pool
from pathlib import Path
from collections import defaultdict
import pysam
from tqdm import tqdm
from field import Variable
from file_io import load_yaml
from run_external import seqtk_subseq


def check_mapped_read(read):
    """Check read mapping flags."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    return True


def calculate_coverage(aln, reads_mapped):
    """Calculate base and read coverage."""
    _base_cov_dict = defaultdict(list)
    read_cov_dict = defaultdict(lambda: 0)
    allowed_operations = set([0, 7, 8])
    with tqdm(total=reads_mapped, unit_scale=True) as pbar:
        for read in aln.fetch(until_eof=True):
            if not check_mapped_read(read):
                continue
            for operation, length in read.cigartuples:
                if operation in allowed_operations:
                    _base_cov_dict[read.reference_name].append(length)
            read_cov_dict[read.reference_name] += 1
            pbar.update()
    base_cov_dict = {ref_name: sum(_base_cov) for ref_name, _base_cov in _base_cov_dict.items()}
    return base_cov_dict, read_cov_dict


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
    stats = {}
    print("Loading mapping data from %s" % bam_file)
    with pysam.AlignmentFile(bam_file, "r%s" % f_char) as aln:
        stats = {'mapped': aln.mapped,
                 'unmapped': aln.unmapped}
        _covs, _read_covs = calculate_coverage(aln, aln.mapped)
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
    fields = {'cov_id': field_id}
    fields['cov_range'] = [min(covs+[kwargs['cov_range'][0]]),
                           max(covs+[kwargs['cov_range'][1]])]
    fields['cov'] = Variable(field_id,
                             values=covs,
                             meta={'field_id': field_id, 'file': bam_file},
                             parents=['children',
                                      {'id': 'base_coverage',
                                       'clamp': 0.01 if fields['cov_range'][0] == 0 else False,
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


def parse_json_cov(json_file, **kwargs):
    """Parse coverage from JSON cov file."""
    data = load_yaml(json_file)
    base_name = Path(json_file).stem.split('.')[0]
    covs = []
    if 'values' in data:
        for value in data['values']:
            covs.append(float("%.4f" % value))
    if base_name.endswith('_read_cov'):
        type = 'read_cov'
        parent = 'read_coverage'
        datatype = 'float'
        clamp = 1
    elif base_name.endswith('_cov'):
        type = 'cov'
        parent = 'base_coverage'
        datatype = 'integer'
        clamp = 0.01
    else:
        return None
    field_id = base_name
    fields = {}
    fields["%s_id" % type] = field_id
    fields["%s_range" % type] = [min(covs+[kwargs["%s_range" % type][0]]),
                                 max(covs+[kwargs["%s_range" % type][1]])]
    if kwargs['meta'].has_field(field_id):
        file_name = kwargs['meta'].field_meta(field_id)['file']
    else:
        file_name = json_file
    fields[type] = Variable(field_id,
                            values=covs,
                            meta={'field_id': field_id, 'file': file_name},
                            parents=['children',
                                     {'id': parent,
                                      'datatype': 'integer',
                                      'clamp': clamp if fields["%s_range" % type][0] == 0 else False,
                                      'range': fields["%s_range" % type]},
                                     'children']
                            )
    return fields


def parse(files, **kwargs):
    """Parse all BAM files."""
    parsed = []
    if kwargs['meta'].has_field('base_coverage'):
        cov_range = kwargs['meta'].field_meta('base_coverage')['range']
    else:
        cov_range = [math.inf, -math.inf]
    if kwargs['meta'].has_field('read_coverage'):
        read_cov_range = kwargs['meta'].field_meta('read_coverage')['range']
    else:
        read_cov_range = [math.inf, -math.inf]
    for file in files:
        if file.endswith('.json'):
            fields = parse_json_cov(file, **kwargs, cov_range=cov_range, read_cov_range=read_cov_range)
        else:
            fields = parse_bam(file, **kwargs, cov_range=cov_range, read_cov_range=read_cov_range)
        if 'cov' in fields:
            parsed.append(fields['cov'])
            cov_range = fields['cov_range']
            if 'y' not in kwargs['meta'].plot:
                kwargs['meta'].plot.update({'y': fields['cov_id']})
        if 'read_cov' in fields:
            parsed.append(fields['read_cov'])
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
