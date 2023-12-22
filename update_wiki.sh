#!/bin/bash

CURRENT_VERSION=$(grep current_version .bumpversion.cfg | head -n 1 | cut -d' ' -f 3)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
WIKI=$SCRIPT_DIR/../blobtoolkit.wiki

update_block () {
    MARKER=$1
    shift
    FILE=$1
    FILE_PATH=$WIKI/$FILE
    shift
    awk -v FLAGS="$@" '/^```'"$MARKER"'/{del=1;print;system("blobtools "FLAGS);next} {if(!del)print} /^```/{if(del){del=0;print}}' \
      $FILE_PATH > $FILE_PATH.tmp && mv $FILE_PATH.tmp $FILE_PATH
}

update_block "sh help text" blobtools-create.md "create -h"
update_block "sh help text" blobtools-add.md "add -h"
update_block "sh help text" blobtools-filter.md "filter -h"
update_block "sh help text" blobtools-host.md "host -h"
update_block "sh help text" blobtools-view.md "view -h"

# sed -E "s:download/[[:digit:]]+\.[[:digit:]]+\.[[:digit:]]+/blobtk-linux:download/v$CURRENT_VERSION/blobtk-linux:; \
#      s:download/[[:digit:]]+\.[[:digit:]]+\.[[:digit:]]+/blobtk-macos-x64:download/v$CURRENT_VERSION/blobtk-macos-x64:; \
#      s:download/[[:digit:]]+\.[[:digit:]]+\.[[:digit:]]+/blobtk-macos-amd64:download/v$CURRENT_VERSION/blobtk-macos-amd64:" \
#      $WIKI/Home.md > $WIKI/Home.md.tmp

