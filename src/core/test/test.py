#!/usr/bin/env python3

from blobtoolkit_core import filter

read_count = filter.fastx(
    {
        "list": ["DWSF010000006.1", "DWSF010000016.1"],
        "bam": "test/DWSF01.SRR12171171.bam",
        "read_list": "test/dict.txt",
    }
)

print(read_count)

read_count = filter.fastx(
    {
        "list_file": "test/test.list",
        "bam": "test/test.bam",
        "fastq1": "test/reads_1.fq.gz",
        "fastq2": "test/reads_2.fq.gz",
        "fastq_out": True,
    }
)

print(read_count)

options = filter.FilterOptions(
    None,
    "test/alt.list",
    "test/DWSF01.SRR12171171.bam",
    None,
    "test/GCA_018252495.1.fasta.gz",
    "test/SRR12171171.fastq.gz",
    None,
    "subset",
    False,
    True,
)
read_count = filter.fastx_with_options(options)

print(read_count)
