#!/bin/bash

# Script to generate BTK plots
# Usage: BLOBDIR=/path/to/blobdir ./generate_btk_plots.sh  

# Path to temporary file directory
if [ -z "$TMPDIR" ]; then
  TMPDIR=/tmp
fi

# Conda environment
if [ -z "$CONDA_ENV" ]; then
  CONDA_ENV=btk_env
fi

# Function to run BTK Viewer as a background process
run_viewer () {
    # Usage: run_viewer <BLOBDIR> <TMPFILE>
    local BLOBDIR=$1
    local TMPFILE=$2
    PYTHONUNBUFFERED=true blobtools view --local $BLOBDIR &>$TMPFILE &
    echo "$!"
}

# Function to list descendant processes
list_descendants () {
  # Usage: stop_viewer <PID> <TMPFILE>
  local PARENT_PID=$1
  local CHILDREN=$(pgrep -P $PARENT_PID)
  for CHILD_PID in $CHILDREN; do
    list_descendants $CHILD_PID
  done
  if [ ! -z "$CHILDREN" ]; then
    echo "$CHILDREN"
  fi
}

# Function to stop BTK Viewer background process
stop_viewer () {
    # Usage: stop_viewer <PID> <TMPFILE>
    local PID=$1
    local TMPFILE=$2
    for CHILD_ID in $(list_descendants $PID); do
        kill -INT $CHILD_ID
    done
    kill -INT $PID
    echo "BTK viewer output:"
    cat $TMPFILE
    rm -f $TMPFILE
}

# Check a valid BlobDir directory has been specified
if [ -z "$BLOBDIR" ]; then
  echo "ERROR: Environment variable BLOBDIR has not been set"
  exit 1
fi
if [ "$BLOBDIR" == "_" ]; then
  # Use example dataset
  DATASET=FXWY01
else
  if [ ! -d "$BLOBDIR" ]; then
    echo "ERROR: $BLOBDIR is not a valid directory path"
    exit 1
  fi
  if [ ! -f "$BLOBDIR/meta.json" ]; then
    echo "ERROR: $BLOBDIR is not a valid BlobDir"
    exit 1
  fi
  # Set dataset name from BLOBDIR
  DATASET=$(basename $BLOBDIR)
fi

# Load the Conda environment
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV
EXITCODE=$?
if [ ! $EXITCODE -eq 0 ]; then
  echo "ERROR: Unable to activate conda environment '$CONDA_ENV'" >&2
  exit 1
fi

# Create temporary file using process ID of this script
TMPFILE="$TMPDIR/$$.out"

# Start Viewer and capture subshell process ID
echo "Starting BTK Viewer"
PID=$(run_viewer $BLOBDIR $TMPFILE)
echo "Waiting for BTK Viewer to start with PID=$PID"

# Parse TMPFILE to find BTK Viewer instance HOST
i=0
while true; do
    HOST=$(grep -oe "http://localhost:80[0-9][0-9]" ${TMPFILE})
    if [ ! -z "$HOST" ]; then
      break
    fi
    if [ $i -gt 9 ]; then
      echo "ERROR: Unable to start viewer after 10s"
      stop_viewer $PID $TMPFILE
      exit 1
    fi
    sleep 1
    i=$((i+1))
done

# Generate BTK plots
echo "Generating plots using BTK Viewer instance at $HOST"

blobtools view --host $HOST $DATASET
blobtools view --host $HOST --view snail $DATASET
blobtools view --host $HOST --view cumulative $DATASET

# Shut down Viewer
echo "Finished generating plots"
echo "Shutting down BTK Viewer instance at $HOST"
stop_viewer $PID $TMPFILE