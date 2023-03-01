"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run BlobToolKit sub-pipelines
-----------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s blobtoolkit.smk
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
  "%s/blastn.stats" % parent_dir,
  "%s/busco.stats" % parent_dir,
  "%s/minimap.stats" % parent_dir,
  "%s/windowmasker.stats" % parent_dir,
  "%s/diamond.stats" % parent_dir,
  "%s/diamond_blastp.stats" % parent_dir,
  "%s/blobtools.stats" % parent_dir,
  "%s/view.stats" % parent_dir
]

# generate_static = config.get("generate_static", False)
# if generate_static:
#   inputs.append("%s/view.stats" % parent_dir)

rule all:
    """
    Dummy rule to define output
    """
    input: inputs

include: "rules/run_sub_pipeline.smk"
