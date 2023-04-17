#!/bin/bash

ls /tmp/*.whl | while read WHEEL; do
    $CONDA_DIR/envs/btk_env/bin/pip install $WHEEL
done

if [ -s /blobtoolkit/blobtoolkit-api-linux ]; then
    mv /blobtoolkit/blobtoolkit-api-linux /blobtoolkit/blobtoolkit-api
    chmod 755 /blobtoolkit/blobtoolkit-api
fi

if [ -s /blobtoolkit/blobtoolkit-viewer-linux ]; then
    mv /blobtoolkit/blobtoolkit-viewer-linux /blobtoolkit/blobtoolkit-viewer
    chmod 755 /blobtoolkit/blobtoolkit-viewer
fi

chown blobtoolkit:blobtoolkit /blobtoolkit/blobtoolkit-{api,viewer}