#!/bin/bash

LEVEL=$1

if [ -z "$LEVEL" ]; then
  echo "Usage: ./bump_version.sh major|minor|patch|prepatch"
  exit 1;
fi

if output=$(git status --porcelain) && [ ! -z "$output" ]; then
  # Uncommitted changes
  echo "Unable to bump version. Git working directory is not clean."
  exit 1
fi

LATEST_TAG=$(grep current_version .prebumpversion.cfg | cut -d' ' -f 3)
NPM_UPDATE=0
HOST_UPDATE=0
PIPELINE_UPDATE=0
BLOBTOOLS_UPDATE=0

if [[ "$LEVEL" == *major ]] || [[ "$LEVEL" == *minor ]]; then
  NPM_UPDATE=1
  HOST_UPDATE=1
  PIPELINE_UPDATE=1
  BLOBTOOLS_UPDATE=1
else
  if [ "$NPM_UPDATE" == 0 ]; then
    git diff --exit-code -s ${LATEST_TAG}..HEAD src/api
    NPM_UPDATE=$?
  fi
  if [ "$NPM_UPDATE" == 0 ]; then
    git diff --exit-code -s ${LATEST_TAG}..HEAD src/viewer
    NPM_UPDATE=$?
  fi
  if [ "$NPM_UPDATE" == 0 ]; then
    git diff --exit-code -s ${LATEST_TAG}..HEAD src/packaged-viewer
    NPM_UPDATE=$?
  fi

  if [ "$HOST_UPDATE" == 0 ]; then
    git diff --exit-code -s ${LATEST_TAG}..HEAD src/blobtoolkit-host
    HOST_UPDATE=$?
    if [ "$HOST_UPDATE" == 1 ]; then
      BLOBTOOLS_UPDATE=1
    fi
  fi

  if [ "$PIPELINE_UPDATE" == 0 ]; then
    git diff --exit-code -s ${LATEST_TAG}..HEAD src/blobtoolkit-pipeline
    PIPELINE_UPDATE=$?
    if [ "$PIPELINE_UPDATE" == 1 ]; then
      BLOBTOOLS_UPDATE=1
    fi
  fi

  if [ "$BLOBTOOLS_UPDATE" == 0 ]; then
    git diff --exit-code -s ${LATEST_TAG}..HEAD src/btk
    BLOBTOOLS_UPDATE=$?
  fi

  if [ "$BLOBTOOLS_UPDATE" == 0 ]; then
    git diff --exit-code -s ${LATEST_TAG}..HEAD src/blobtools
    BLOBTOOLS_UPDATE=$?
  fi
  
  if [ "$BLOBTOOLS_UPDATE" == 0 ]; then
    LEVEL="-preid rc pre$LEVEL"
  fi
fi

if [ "$BLOBTOOLS_UPDATE" == 1 ]; then
  bump2version --allow-dirty --config-file .prebumpversion.cfg $LEVEL
  VERSION=$(grep current_version .prebumpversion.cfg | cut -d' ' -f 3)
  if [ "$HOST_UPDATE" == 1 ]; then
    bump2version --allow-dirty --config-file .hostbumpversion.cfg $VERSION
  fi
  if [ "$PIPELINE_UPDATE" == 1 ]; then
    bump2version --allow-dirty --config-file .hostbumpversion.cfg $VERSION
  fi
  ./bump_full_version.sh $NPM_UPDATE $VERSION
elif [ "$NPM_UPDATE" == 1 ]; then
  ./bump_npm_version.sh $LEVEL
fi



