---
layout: default
title: Command Line Tutorials
---

# Command Line Tutorials

Step-by-step walkthroughs for common BlobTools command-line tasks.

## Workflow 1: Create a Dataset

The minimum input for dataset creation is an assembly FASTA file. Metadata and taxonomy can be added at creation time.

### Step 1: Prepare Your Data

You'll need:

- Assembly FASTA file
- Optional metadata file (YAML/JSON)
- Optional taxonomy context (`--taxid` with `--taxdump`)

### Step 2: Create BlobDir

```bash
blobtools create --fasta assembly.fasta \
                 --meta metadata.json \
                 my_blobdir/
```

Example including taxonomy:

```bash
blobtools create \
    --fasta assembly.fasta \
    --meta metadata.yaml \
    --taxid 75913 \
    --taxdump /path/to/taxdump \
    my_blobdir/
```

## Workflow 2: Add Analyses

Use `blobtools add` to import additional analyses after creation.

### Step 3: Add Data

```bash
blobtools add --hits blast_results.txt \
              --cov coverage.sorted.bam \
              my_blobdir/
```

You can also import BUSCO and text-based fields:

```bash
blobtools add \
    --busco full_table.tsv \
    --text custom_metrics.tsv \
    my_blobdir/
```

Use `blobtools replace` if existing values should be overwritten.

## Workflow 3: Filter and Export

Filter by parameter values and write a filtered dataset or files.

### Step 4: Filter Using Parameters

```bash
blobtools filter \
    --param length--Min=1000 \
    --param bestsumorder_phylum--Keys=no-hit \
    --summary STDOUT \
    my_blobdir/
```

### Step 5: Reproduce Viewer Filters on CLI

```bash
blobtools filter \
    --query-string "gc--Min=0.3&bestsumorder_phylum--Keys=no-hit" \
    --fasta assembly.fasta \
    my_blobdir/
```

Or load an exported viewer selection/list:

```bash
blobtools filter \
    --json list.json \
    --fasta assembly.fasta \
    my_blobdir/
```

### Step 6: Visualize Locally

```bash
blobtools host --port 8080 my_blobdir/
```

Then visit `http://localhost:8080` in your browser.

## Further Reading

See the [Commands Reference](commands) for complete option details.
