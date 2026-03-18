---
layout: default
title: Installation
---

# Installation

BlobToolKit can be installed using Conda, pip, or Docker.

## Conda Installation

The easiest way to install BlobToolKit is using Conda:

```bash
conda install -c bioconda blobtoolkit
```

## Pip Installation

You can also install BlobToolKit using pip:

```bash
pip install blobtoolkit
```

To include optional command groups:

```bash
pip install "blobtoolkit[host]"
pip install "blobtoolkit[pipeline]"
pip install "blobtoolkit[full]"
```

## Docker Installation

A Docker container is available for containerized deployments:

```bash
docker pull genomehubs/blobtoolkit:latest
docker run -it genomehubs/blobtoolkit:latest blobtools --version
```

## Optional Viewer Dependencies

For some host/view workflows you may need browser and display tooling (for example Firefox and X11/Xvfb depending on platform and runtime mode).

## Verify Installation

To verify your installation works correctly:

```bash
blobtools --version
btk --version
```

## System Requirements

- Python 3.10 or higher
- 2GB RAM minimum
- 10GB free disk space (for databases)

## Databases

For taxonomy and similarity workflows, prepare local database resources:

- NCBI taxdump
- NCBI nt
- UniProt reference proteomes
- BUSCO lineage data

See [Database Setup](databases) for commands and structure.

## Next Steps

After installation, check out the [Pipeline Getting Started](../pipeline/getting-started) guide.
