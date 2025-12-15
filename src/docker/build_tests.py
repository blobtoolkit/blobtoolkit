#!/usr/bin/env python3

"""Run a set of simple tests to check if the commands are available and working."""

import subprocess
import sys

commands = [
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
    ["btk", "pipeline", "run", "--help"],
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
        subprocess.run(command, check=True)
        print(f"Command {' '.join(command)}: SUCCESS")
    except subprocess.CalledProcessError:
        print(f"Command {' '.join(command)}: FAIL")
        return False
    return True


def main():
    all_success = True
    for command in commands:
        if not run_command(command):
            all_success = False

    if not all_success:
        sys.exit(1)


if __name__ == "__main__":
    main()
