---
layout: default
title: Command Line Tutorials
---

# Command Line Tutorials

Step-by-step walkthroughs for common command line tasks.

## Creating a Basic Blob Plot

This tutorial walks through creating a simple blob plot dataset.

### Step 1: Prepare Your Data

You'll need:

- Assembly FASTA file
- Taxonomic assignment results (BLAST/Diamond)
- Coverage information (BAM files)

### Step 2: Create BlobDir

```bash
blobtools create --fasta assembly.fasta \
                 --meta metadata.json \
                 my_blobdir/
```

### Step 3: Add Data

```bash
blobtools add --hits blast_results.txt \
              --cov coverage.sorted.bam \
              my_blobdir/
```

### Step 4: Visualize

```bash
blobtools host --port 8080 my_blobdir/
```

Then visit `http://localhost:8080` in your browser.

## Filtering Contigs

Use the filter command to subset your data:

```bash
blobtools filter --query "read_cov > 5" my_blobdir/ --output_dir filtered/
```

## Further Reading

See the [Commands Reference](commands.md) for more options. Visit the [Wiki](https://github.com/genomehubs/blobtoolkit/wiki) for additional tutorials.
