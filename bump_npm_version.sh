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

# Update Dockerfiles
VERSION=$(grep '"version":' src/api/package.json | cut -d'"' -f 4)
sed -i '/ARG VERSION=/s/.*/ARG VERSION='$VERSION'/' src/docker/api/Dockerfile
sed -i '/ARG VERSION=/s/.*/ARG VERSION='$VERSION'/' src/docker/viewer/Dockerfile

git commit -a -m "bump UI/API version"
