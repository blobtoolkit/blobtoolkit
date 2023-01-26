import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run BUSCO
--------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s busco.smk
    -j 8

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2022 Genome Research Limited, MIT License
"""

include: "../lib/functions.py"

chunk_stats_path = "../chunk_stats"

rule all:
    """
    Dummy rule to define all outputs
    """
    input:
        expand("%s.busco.{lineage}/full_table.tsv.gz" % config["assembly"]["prefix"], lineage=config["busco"]["lineages"]),
        "%s.chunk_stats.tsv"  % config["assembly"]["prefix"],


include: "rules/run_busco5.smk"
include: "rules/count_busco_genes.smk"
include: "rules/unzip_assembly_fasta.smk"
