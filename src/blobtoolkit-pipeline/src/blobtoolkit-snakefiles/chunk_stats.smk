import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to generate sequence stats for sequence chunks
-------------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s chunk_stats.smk
    -j 8

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2022 Genome Research Limited  % config["assembly"]["prefix"], MIT License
"""

include: "../lib/functions.py"

minimap_path = "../minimap"
windowmasker_path = "../windowmasker"

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s.chunk_stats.mask.bed"  % config["assembly"]["prefix"],
        "%s.chunk_stats.tsv"  % config["assembly"]["prefix"],

include: "rules/get_chunked_stats.smk"
