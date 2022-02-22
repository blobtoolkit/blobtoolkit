#!/bin/bash

cd ./src/packaged-viewer &&

./package.sh &&

cd - &&

mkdir -p dist &&

rm -rf dist/* &&

mv ./src/packaged-viewer/dist/blobtoolkit-viewer* ./dist/ &&

mv ./src/viewer/dist/blobtoolkit-viewer* ./dist/
