#!/bin/bash

# Script to generate BTK plots
# Usage: BLOBDIR=/path/to/blobdir ./generate_btk_plots.sh

# Path to temporary file directory
if [ -z "$TMPDIR" ]; then
  TMPDIR=/tmp
fi
TMPDIR="$TMPDIR/btk-$$"

# Conda environment
if [ -z "$CONDA_ENV" ]; then
  CONDA_ENV=btk_env
fi

# Timeout to wait for farm latency
if [ -z "$TIMEOUT" ]; then
  TIMEOUT=180
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
  # Usage: stop_viewer <PID> <TMPFILE> <TMPDIR> <MAINDIR>
  local PID=$1
  local TMPFILE=$2
  local TMPDIR=$3
  local MAINDIR=$4
  echo "Stopping child processes for $PID"
  for CHILD_ID in $(list_descendants $PID); do
    echo " - stopping $CHILD_ID"
    kill -INT $CHILD_ID
    while kill -0 $PID; do 
      sleep 2
      kill -0 $PID && kill -9 $PID
    done
    echo " - stopped $CHILD_ID"
  done
  echo " - stopping $PID"
  kill -INT $PID
  while kill -0 $PID; do
    sleep 2
    kill -0 $PID && kill -9 $PID
  done
  echo " - all processes stopped"
  if [ ! -z "$MAINDIR" ]; then
    echo "Moving files from $TMPDIR to $MAINDIR"
    mv "$TMPDIR/*.png" "$MAINDIR/"
  else
    echo "BTK viewer output:"
    cat $TMPFILE
  fi
  echo "Removing $TMPDIR"
  rm -rf $TMPDIR
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
  # Set dataset name from BLOBDIR and copy to TMPDIR
  DATASET=$(basename $BLOBDIR)
  mkdir -p "$TMPDIR"
  cp -r $BLOBDIR "$TMPDIR/$DATASET"
  BLOBDIR="$TMPDIR/$DATASET"
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
mkdir -p "$TMPDIR"
TMPFILE="$TMPDIR/$$.out"

# Start Viewer and capture subshell process ID
echo "Starting BTK Viewer for $BLOBDIR"
PID=$(run_viewer $BLOBDIR $TMPFILE)
echo "Waiting for BTK Viewer to start with PID=$PID"

# Parse TMPFILE to find BTK Viewer instance HOST
i=0
while true; do
  i=$((i+1))
  HOST=$(grep -oe "http://localhost:80[0-9][0-9]" ${TMPFILE})
  if [ ! -z "$HOST" ]; then
    break
  fi
  if [ $i -gt $TIMEOUT ]; then
    echo "ERROR: Unable to start viewer after ${TIMEOUT}s"
    stop_viewer $PID $TMPFILE $TMPDIR
    exit 1
  fi
  sleep 1
done

# Wait for BTK Viewer to be available
echo "Checking connection to BTK Viewer instance"
i=0
while true; do
  i=$((i+1))
  if curl --output /dev/null --silent --head --fail "$HOST/view/"; then
    echo "Connected to BTK Viewer instance at $HOST"
    break
  fi
  if [ $i -gt $TIMEOUT ]; then
    echo "ERROR: Unable to connect to viewer after ${TIMEOUT}s"
    stop_viewer $PID $TMPFILE $TMPDIR
    exit 1
  fi
  sleep 1
done

# Generate BTK plots
echo "Generating plots using BTK Viewer instance at $HOST"
MAINDIR=$(pwd)
cd $TMPDIR
PYTHONUNBUFFERED=true timeout $TIMEOUT blobtools view --host $HOST --view blob --param plotShape=circle $DATASET
PYTHONUNBUFFERED=true timeout $TIMEOUT blobtools view --host $HOST --view snail $DATASET
PYTHONUNBUFFERED=true timeout $TIMEOUT blobtools view --host $HOST --view cumulative $DATASET
cd $MAINDIR

# Shut down Viewer
echo "Finished generating plots"
echo "Shutting down BTK Viewer instance at $HOST"
stop_viewer $PID $TMPFILE $TMPDIR $MAINDIR

echo "Finished script"
exit 0