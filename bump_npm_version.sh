#!/bin/bash

LEVEL=$1

if [ -z "$LEVEL" ]; then
  echo "Usage: ./bump_version.sh major|minor|patch"
  exit 1;
fi

echo $LEVEL
exit

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

git commit -a -m "bump UI/API version"
