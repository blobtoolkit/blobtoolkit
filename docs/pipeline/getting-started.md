# Pipeline Getting Started

BlobToolKit provides pipeline implementations using Nextflow and Snakemake for automated genome assembly quality evaluation.

## Available Workflows

### Nextflow

The Nextflow pipeline offers:

- Reproducible workflows across different systems
- Scalability to high-performance computing environments
- Integrated error handling and resource management

```bash
nextflow run genomehubs/blobtoolkit-pipeline
```

### Snakemake

The Snakemake pipeline provides:

- Flexible Python-based workflow definition
- Integration with conda environments
- Fine-grained control over individual steps

```bash
snakemake --snakefile blobtoolkit.smk
```

## Basic Workflow Steps

1. Prepare your assembly and reads
2. Configure the pipeline (see [Configuration](configuration.md))
3. Run the pipeline
4. View results in the [BlobToolKit Viewer](../viewer/local-hosting.md)

## Next Steps

- Learn about pipeline configuration in [Configuration Guide](configuration.md)
- Check the [command line](../command-line/commands.md) reference
