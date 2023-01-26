import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to generate sequence stats across various window sizes
---------------------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s window_stats.smk
    -j 8

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2022 Genome Research Limited  % config["assembly"]["prefix"], MIT License
"""

include: "../lib/functions.py"

cov_stats_path = "../cov_stats"

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s.window_stats.tsv" % config["assembly"]["prefix"],

include: "rules/get_window_stats.smk"
