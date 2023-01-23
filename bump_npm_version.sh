#!/bin/bash

LEVEL=$1

if [ -z "$LEVEL" ]; then
  echo "Usage: ./bump_version.sh major|minor|patch"
  exit 1;
fi

CURRENT_VERSION=$(grep '"version":' src/viewer/package.json | cut -d'"' -f 4)

if [[ "$LEVEL" == pre* ]]; then
  LEVEL="-preid rc $LEVEL"
fi

cd src/api &&

npm version --no-git-tag-version $LEVEL &&

cd - &&

cd src/viewer &&

npm version --no-git-tag-version $LEVEL &&

cd - &&

cd src/packaged-viewer &&

npm version --no-git-tag-version $LEVEL &&

cd -

NEW_VERSION=$(grep '"version":' src/viewer/package.json | cut -d'"' -f 4)

git commit -a -m "bump UI/API version: ${CURRENT_VERSION} → ${NEW_VERSION}"
git tag -a $NEW_VERSION -m "Bump version: ${CURRENT_VERSION} → ${NEW_VERSION}"
