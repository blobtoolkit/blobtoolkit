#!/bin/bash

LEVEL=$1
VERSION=$2
NPM_UPDATE=$3

CURRENT_VERSION=$(grep current_version .bumpversion.cfg | head -n 1 | cut -d' ' -f 3)

if [ $NPM_UPDATE == 1 ]; then
  cd src/api &&
  npm version --no-git-tag-version $VERSION &&
  cd - &&

  cd src/viewer &&
  npm version --no-git-tag-version $VERSION &&
  cd - &&

  cd src/packaged-viewer &&
  npm version --no-git-tag-version $VERSION &&
  cd -
fi

bump2version --allow-dirty $LEVEL

NEW_VERSION=$(grep current_version .bumpversion.cfg | head -n 1 | cut -d' ' -f 3)

git commit -a -m "Bump version: ${CURRENT_VERSION} → ${NEW_VERSION}"
git tag -a $NEW_VERSION -m "Bump version: ${CURRENT_VERSION} → ${NEW_VERSION}"
