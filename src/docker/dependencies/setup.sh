#!/bin/bash

ls /tmp/*.whl | while read WHEEL; do
    echo y | $CONDA_DIR/envs/btk_env/bin/pip uninstall ${WHEEL%%-*}
    $CONDA_DIR/envs/btk_env/bin/pip install $WHEEL
done

# touch /tmp/build_env.t$BUILD_ENV

# # if [ "${BUILD_ENV}" == "dev" ]; then
#   cd /tmp/blobtoolkit
#   git fetch origin feature/split-packaging
#   git checkout feature/split-packaging

#   pip install setuptools wheel twine

#   cd /tmp/blobtoolkit/src/blobtoolkit-pipeline
#   ./pip_install_latest.sh manylinux2014_x86_64

#   cd /tmp/blobtoolkit/src/blobtoolkit-host
#   ./pip_install_latest.sh manylinux2014_x86_64

#   cd /tmp/blobtoolkit
#   ./pip_install_latest.sh manylinux2014_x86_64
# # fi