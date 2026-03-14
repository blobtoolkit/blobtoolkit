#!/bin/bash

ls /tmp/*.whl | while read WHEEL; do
    $CONDA_DIR/envs/btk_env/bin/pip install \
        --prefix $CONDA_DIR/envs/btk_env \
        --force-reinstall \
        $WHEEL
done

# Fall back to PyPI for host/pipeline if no local wheel was built for them
$CONDA_DIR/envs/btk_env/bin/python -c "import importlib.metadata; importlib.metadata.version('blobtoolkit-host')" 2>/dev/null || \
    $CONDA_DIR/envs/btk_env/bin/pip install --prefix $CONDA_DIR/envs/btk_env "blobtoolkit[host]"

$CONDA_DIR/envs/btk_env/bin/python -c "import importlib.metadata; importlib.metadata.version('blobtoolkit-pipeline')" 2>/dev/null || \
    $CONDA_DIR/envs/btk_env/bin/pip install --prefix $CONDA_DIR/envs/btk_env "blobtoolkit[pipeline]"

if [ -s /blobtoolkit/blobtoolkit-api-linux ]; then
    mv /blobtoolkit/blobtoolkit-api-linux /blobtoolkit/blobtoolkit-api
    chmod 755 /blobtoolkit/blobtoolkit-api
fi

if [ -s /blobtoolkit/blobtoolkit-viewer-linux ]; then
    mv /blobtoolkit/blobtoolkit-viewer-linux /blobtoolkit/blobtoolkit-viewer
    chmod 755 /blobtoolkit/blobtoolkit-viewer
fi

chown blobtoolkit:blobtoolkit /blobtoolkit/blobtoolkit-{api,viewer}