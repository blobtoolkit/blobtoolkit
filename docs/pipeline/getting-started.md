---
layout: default
title: Pipeline Getting Started
---

# Pipeline Getting Started

The preferred pipeline for BlobToolKit is the Nextflow workflow documented at:

- https://pipelines.tol.sanger.ac.uk/blobtoolkit

## Preferred Workflow

Use the external pipeline docs for:

- configuration format and parameters
- execution instructions
- platform-specific guidance
- current recommended defaults

At present, this preferred Nextflow workflow is configured and run outside the BlobToolKit CLI.

Integration of this preferred workflow into `btk pipeline run` may be added in a future update.

## Legacy/Local CLI Workflow

The local CLI pipeline entry points are still available for internal or legacy usage:

```bash
btk pipeline data --help
btk pipeline run --help
```

These commands dispatch to the installed blobtoolkit-pipeline package.

## Alternative Workflows

### Nextflow

Nextflow remains the preferred approach. Use:

- https://pipelines.tol.sanger.ac.uk/blobtoolkit

rather than older local run instructions.

### Snakemake (legacy/advanced)

Some historical workflows were run directly through Snakemake. This can still be useful in advanced HPC environments, but should be treated as an alternative path when maintaining older setups.

```bash
snakemake -p --use-conda --configfile /path/to/assembly.yaml
```

## Basic Workflow Steps

1. Prepare your assembly and reads
2. Follow configuration instructions at https://pipelines.tol.sanger.ac.uk/blobtoolkit
3. Run the preferred Nextflow workflow as documented there
4. View results in the [BlobToolKit Viewer](../viewer/local-hosting)

## Next Steps

- Review local/legacy config fields in [Configuration Guide](configuration)
- Check the [command line](../command-line/commands) reference
