#!/bin/bash

echo "Packaging viewer"

TEMPLATE="<% if (variables) { %> \
<script> \
<%- variables %> \
</script> \
<% } %> \
"

cp .env.dist .env &&

npm run build &&

rm -rf ./ui/src/public &&

tar -C ./dist -czf ./dist/blobtoolkit-viewer.tgz public &&

mv ./dist/public ./ui/src/ &&

mkdir -p ./ui/src/views &&

rm -rf ./ui/src/views/* &&

sed 's:<!---->:'"$TEMPLATE"':' ./ui/src/public/index.html > ./ui/src/views/index.ejs &&

pkg --compress GZip ui/package.json &&

cp ./api/.env.dist ./api/.env &&

pkg --compress GZip api/package.json