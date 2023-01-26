import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run Windowmasker
--------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s windowmasker.smk
    -j 8

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2022 Genome Research Limited, MIT License
"""

include: "lib/functions.py"

rule all:
    """
    Dummy rule to define all outputs
    """
    input:
        "%s.windowmasker.fasta" % config["assembly"]["prefix"]


include: "rules/run_windowmasker.smk"
include: "rules/unzip_assembly_fasta.smk"
