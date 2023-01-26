import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run blastn
--------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s blastn.smk
    -j 8

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2022 Genome Research Limited, MIT License
"""

include: "../lib/functions.py"

diamond_path = "../diamond"
windowmasker_path = "../windowmasker"

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s.blastn.nt.out" % config["assembly"]["prefix"]

include: "rules/extract_nohit_fasta.smk"
include: "rules/chunk_nohit_fasta.smk"
include: "rules/run_blastn.smk"
include: "rules/unchunk_blastn.smk"
