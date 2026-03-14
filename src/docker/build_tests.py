#!/usr/bin/env python3

"""Run a set of simple tests to check if the commands are available and working."""

import subprocess
import sys
from importlib import metadata

required_commands = [
    ["blobtools", "--version"],
    ["blobtools", "--help"],
    ["blobtools", "create", "--help"],
    ["blobtools", "replace", "--help"],
    ["blobtools", "add", "--help"],
    ["blobtools", "remove", "--help"],
    ["blobtools", "validate", "--help"],
    ["blobtools", "filter", "--help"],
    ["blobtools", "host", "--help"],
    ["blobtools", "view", "--help"],
    ["btk", "pipeline", "--help"],
    ["btk", "pipeline", "add-summary-to-metadata", "--help"],
    ["btk", "pipeline", "chunk-fasta", "--help"],
    ["btk", "pipeline", "count-busco-genes", "--help"],
    ["btk", "pipeline", "extract-busco-genes", "--help"],
    ["btk", "pipeline", "generate-config", "--help"],
    ["btk", "pipeline", "generate-static-images", "--help"],
    ["btk", "pipeline", "transfer-completed", "--help"],
    ["btk", "pipeline", "unchunk-blast", "--help"],
    ["btk", "pipeline", "window-stats", "--help"],
]


def run_command(command):
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"Command {' '.join(command)}: SUCCESS")
    except subprocess.CalledProcessError as exc:
        print(f"Command {' '.join(command)}: FAIL")
        if exc.stdout:
            print(exc.stdout.strip())
        if exc.stderr:
            print(exc.stderr.strip())
        return False
    return True


def require_installed(package):
    try:
        version = metadata.version(package)
        print(f"Package {package}=={version}: INSTALLED")
        return True
    except metadata.PackageNotFoundError:
        print(f"Package {package}: MISSING")
        return False


def main():
    all_success = True

    if not require_installed("blobtoolkit-host"):
        all_success = False
    if not require_installed("blobtoolkit-pipeline"):
        all_success = False

    for command in required_commands:
        if not run_command(command):
            all_success = False

    if not all_success:
        sys.exit(1)


if __name__ == "__main__":
    main()
