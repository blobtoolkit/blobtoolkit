#!/bin/bash

if [[ ! -z $VIEWER ]]; then
  blobtools view /blobtoolkit/datasets &
  tail -f /dev/null
  exit
fi

if [[ ! -z $ASSEMBLY ]]; then

  export WORKDIR=/blobtoolkit/datasets/$ASSEMBLY

  if [[ ! -f "$WORKDIR/config.yaml" ]]; then
    echo "ERROR: unable to locate configuration file '$WORKDIR/config.yaml'"
    exit 1
  fi
  
  if [[ -z $THREADS ]]; then
    export THREADS=32
  fi

  if [[ ! -z $DRYRUN ]]; then
    DRYRUN="--dry-run"
  fi

  if [[ -z $TOOL ]]; then
    TOOL="blobtoolkit"
  fi

  if [[ ! -z $UNLOCK ]]; then
    btk pipeline run --unlock --config $WORKDIR/config.yaml
  fi

  btk pipeline run $DRYRUN --threads $THREADS --config $WORKDIR/config.yaml

  if [ $? -ne 0 ];then
    echo "ERROR: failed to run snakemake pipeline"
    exit 1
  fi

  if [[ ! -z $TRANSFER ]]; then
    btk pipeline transfer-completed \
        --in $WORKDIR \
        --out /blobtoolkit/output
    
    if [ $? -ne 0 ];then
      echo "ERROR: failed to transfer completed files"
      exit 1
    fi
  fi
fi


