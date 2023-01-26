"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run add reads to an existing BTK run
------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s add_reads.smk
    -j 60

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2022 Genome Research Limited, MIT License
"""

import os

include: '../lib/functions.py'

working_dir = os.getcwd()
parent_dir = os.path.dirname(working_dir)

inputs = [
  "%s/minimap.stats" % parent_dir,
  "%s/blobtools.stats" % parent_dir,
  "%s/view.stats" % parent_dir
]

rule all:
    """
    Dummy rule to define output
    """
    input: inputs

include: "rules/run_sub_pipeline.smk"
