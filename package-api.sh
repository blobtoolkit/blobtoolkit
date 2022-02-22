#!/bin/bash

cd ./src/api &&

./package.sh &&

cd - &&

mkdir -p dist &&

mv ./src/api/dist/genomehubs-api* ./dist/