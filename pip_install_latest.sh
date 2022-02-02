#!/bin/bash

# Build and install the latest blobtools2 version in the current environment

BTK_VERSION=$(
    grep current_version `dirname "$0"`/.bumpversion.cfg \
    | head -n 1 \
    | awk '{print $3}')

python3 setup.py sdist bdist_wheel &&
ls dist &&
echo y | pip uninstall blobtools2 &&
pip install dist/blobtools2-${BTK_VERSION}-py3-none-any.whl ||
echo failed
