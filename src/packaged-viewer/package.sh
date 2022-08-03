#!/bin/bash

echo "Packaging blobtoolkit-viewer" &&

cd ../viewer &&

cp .env.dist .env &&

echo "Installing dependencies" &&

npm install --legacy-peer-deps &&

echo "Bundling javascript" &&

npm run build &&

echo "Archiving build"

tar -C dist -czf dist/blobtoolkit-viewer.tgz public &&

cd - &&

echo "Preparing template" &&

mkdir -p src/views &&

rm -rf src/public &&

mv ../viewer/dist/public src/public &&

rm -rf src/views/* &&

TEMPLATE="<% if (variables) { %> \
<script> \
<%- variables %> \
</script> \
<% } %> \
" &&

sed 's:<!---->:'"$TEMPLATE"':' src/public/index.html > src/views/index.ejs &&

rm src/public/index.html &&

echo "Installing package dependencies" &&

npm install &&

echo "Creating package" &&

pkg --compress GZip package.json
