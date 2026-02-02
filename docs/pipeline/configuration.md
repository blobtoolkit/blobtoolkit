# Pipeline Configuration

Configure the BlobToolKit pipeline using YAML configuration files and samplesheets.

## Configuration File (config.yaml)

The main configuration file controls pipeline parameters:

```yaml
# Genome assembly
assembly: path/to/assembly.fasta

# Reference database
reference: /path/to/reference.fasta

# BUSCO lineage
busco_lineage: lineage_dataset.tar.gz

# Output directory
output_dir: results/

# Threads and memory
threads: 8
memory: 16
```

## Samplesheet Format

Define your samples in a samplesheet file:

```
sample,reads1,reads2,type
sample1,reads1_R1.fq,reads1_R2.fq,paired
sample2,reads2_R1.fq,reads2_R2.fq,paired
sample3,reads3.fq,,single
```

## Running with Configuration

```bash
nextflow run blobtoolkit-pipeline -c config.yaml --samplesheet samplesheet.csv
```

## Configuration Options

Refer to the pipeline documentation for a complete list of configuration options.
