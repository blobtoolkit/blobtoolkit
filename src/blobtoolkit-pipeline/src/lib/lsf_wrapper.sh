#!/bin/bash
#BSUB -n 32
#BSUB -o btk_v2.%J.out
#BSUB -e btk_v2.%J.err
#BSUB -q long
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -M 100000

# Check required variables have been specified
if [ -z "$DATA_DIR" ]; then
  echo "ERROR: you must specify DATA_DIR=/path/to/data" >&2
  exit 1
fi
if [ -z "$SNAKE_DIR" ]; then
  echo "ERROR: you must specify SNAKE_DIR=/path/to/pipeline/directory" >&2
  exit 1
fi
if [ -z "$ACCESSION" ]; then
  echo "ERROR: you must specify ACCESSION, e.g. ACCESSION=GCA_00000000.0" >&2
  exit 1
fi
if [ -z "$TOOL" ]; then
  echo "ERROR: you must specify TOOL, e.g. TOOL=busco" >&2
  exit 1
fi

# Check the configuration file exists
if [ ! -s "$DATA_DIR/$ACCESSION/config.yaml" ]; then
  echo "ERROR: config file $DATA_DIR/$ACCESSION/config.yaml does not exist" >&2
  exit 1
fi

# Check the Snakefile exists
if [ ! -s "$SNAKE_DIR/$TOOL.smk" ]; then
  echo "ERROR: $SNAKE_DIR/$TOOL.smk is not a valid Snakefile" >&2
  exit 1
fi


# Load the Conda environment containing snakemake
eval "$(conda shell.bash hook)"
conda activate bts_env

EXITCODE=$?
if [ ! $EXITCODE -eq 0 ]; then
  echo "ERROR: Unable to activate conda environment 'bts_env'" >&2
  exit 1
fi

# Check the working directory is unlocked in case a previous run failed
mkdir -p $DATA_DIR/$ACCESSION/$TOOL
snakemake -p \
          -j 1 \
          --directory $DATA_DIR/$ACCESSION/$TOOL \
          --configfile $DATA_DIR/$ACCESSION/config.yaml \
          -s $SNAKE_DIR/$TOOL.smk \
          --unlock

if [ $? -ne 0 ];then
  echo "ERROR: failed while unlocking working directory"
  exit 1
fi

# Run pipeline
snakemake -p \
          -j $THREADS \
          --directory $DATA_DIR/$ACCESSION/$TOOL \
          --configfile $DATA_DIR/$ACCESSION/config.yaml \
          --latency-wait 60 \
          --stats $DATA_DIR/$ACCESSION/$TOOL.stats \
          -s $SNAKE_DIR/$TOOL.smk

if [ $? -ne 0 ];then
  echo "ERROR: failed while running pipeline"
  exit 1
fi

echo "done"