#!/bin/bash

# Helper script to generate a configuration file for a public assembly,
# download files and run BTK pipeline

# Check required variables have been specified
ACCESSION=$1
if [ -z "$ACCESSION" ]; then
  echo "ERROR: you must specify ACCESSION, e.g. $0 GCA_00000000.0" >&2
  exit 1
fi

DATA_DIR=/volumes/data/by_accession
DONE_DIR=/volumes/data/processed_by_accession
SNAKE_DIR=~/blobtoolkit/insdc-pipeline
THREADS=60

RESTART=$2

if [ -z "$RESTART" ]; then
  # generate the config file
  $SNAKE_DIR/lib/generate_config.py $ACCESSION \
        --download \
        --out /volumes/data/by_accession \
        --db /volumes/databasenfs \
        --db-suffix 2021_06

  if [ $? -ne 0 ];then
  echo "ERROR: failed to generate config file"
  exit 1
  fi
fi

# Load the Conda environment containing snakemake
eval "$(conda shell.bash hook)"
conda activate bts_env

# Run the snakemake pipeline
mkdir -p $DONE_DIR
TOOL=blobtoolkit
snakemake -p \
          -j $THREADS \
          --directory $DATA_DIR/$ACCESSION/$TOOL \
          --configfile $DATA_DIR/$ACCESSION/config.yaml \
          --latency-wait 60 \
          --stats $DATA_DIR/$ACCESSION/$TOOL.stats \
          -s $SNAKE_DIR/$TOOL.smk

if [ $? -ne 0 ];then
  echo "ERROR: failed to run snakemake pipeline"
  exit 1
fi

if [ ! -z "$RESTART" ]; then

  rm -rf $DATA_DIR/$ACCESSION/reads

  rm -rf $DATA_DIR/$ACCESSION/assembly

  rm -rf $DATA_DIR/$ACCESSION/minimap/*.bam

  mv $DATA_DIR/$ACCESSION $DONE_DIR/$ACCESSION

fi
