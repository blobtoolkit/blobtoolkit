#!/bin/bash

# Build and install the latest blobtoolkit version in the current environment

PLATFORM=$1

if [ -z "$PLATFORM" ]; then
    echo "USAGE: ./pip_install_latest.sh linux_x86_64"
    echo "   OR: ./pip_install_latest.sh macosx_10_9_x86_64"
    echo "   OR: ./pip_install_latest.sh macosx_11_0_arm64"
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# Read version directly from local setup.py so this script works when the
# directory is deployed standalone (e.g. rsync'd without the parent repo).
BTK_VERSION=$(sed -n 's/.*version="\([^"]*\)".*/\1/p' "$SCRIPT_DIR/setup.py" | head -n 1)

if [ -z "$BTK_VERSION" ]; then
    echo "ERROR: Could not determine version from $SCRIPT_DIR/setup.py"
    exit 1
fi

# Use pypa/build if available (avoids setup.py deprecation warnings);
# fall back to legacy setup.py bdist_wheel.
if python3 -m build --version > /dev/null 2>&1; then
    python3 -m build --wheel --no-isolation --config-setting="--plat-name=$PLATFORM"
else
    python3 setup.py bdist_wheel --python-tag py3 --plat-name="$PLATFORM"
fi &&
echo y | pip3 uninstall blobtoolkit-pipeline &&
pip3 install "$SCRIPT_DIR/dist/blobtoolkit_pipeline-${BTK_VERSION}-py3-none-$PLATFORM.whl" &&
blobtoolkit-pipeline -v
