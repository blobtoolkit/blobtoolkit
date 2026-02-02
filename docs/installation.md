---
layout: default
title: Installation
---

# Installation

BlobToolKit can be installed using Conda, Pip, or Docker. Choose the method that best suits your environment.

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

## Docker Installation

A Docker container is available for containerized deployments:

```bash
docker pull genomehubs/blobtoolkit:latest
docker run -it genomehubs/blobtoolkit:latest blobtools --version
```

## Verify Installation

To verify your installation works correctly:

```bash
blobtools --version
```

## System Requirements

- Python 3.8 or higher
- 2GB RAM minimum
- 10GB free disk space (for databases)

## Next Steps

After installation, check out the [Pipeline Getting Started](pipeline/getting-started.md) guide.
