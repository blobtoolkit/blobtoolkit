#!/bin/bash

# Helper script to loop through a file of target accessions

todo="accessions.todo"
active="accessions.active"
failed="accessions.failed"
while read -r line; do
    echo $line
    ./run_btk_pipeline.sh $line < /dev/null &&
    echo $line >> $active ||
    echo $line >> $failed
    sed -i "1 d" $todo
done < $todo
