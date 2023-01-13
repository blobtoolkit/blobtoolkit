#!/bin/bash

NPM_UPDATE=$1

VERSION=$(grep current_version .prebumpversion.cfg | cut -d' ' -f 3)

if [ $NPM_UPDATE == 1 ]; then
  cd src/api &&
  npm version --no-git-tag-version $VERSION &&
  cd - &&

  cd src/viewer &&
  npm version --no-git-tag-version $VERSION &&
  cd - &&

  cd src/packaged-viewer &&
  npm version --no-git-tag-version $VERSION &&
  cd - &&

  git commit -a -m "bump UI/API version"
fi

bump2version $VERSION
