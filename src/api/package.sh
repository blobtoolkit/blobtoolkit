#!/bin/bash

echo "Packaging blobtoolkit-api" &&

echo "Installing dependencies" &&

npm install &&

echo "Copying files" &&

mkdir -p src/api-docs &&

rm -rf src/api-docs/* &&

cp node_modules/swagger-ui-dist/* src/api-docs/ &&

cp node_modules/swagger-ui-express/* src/api-docs/ &&

mkdir -p demo &&

rm -rf demo/* &&

cp -r ../demo/* demo/ &&

cp .env.dist .env &&

echo "Creating package" &&

pkg --compress GZip package.json
