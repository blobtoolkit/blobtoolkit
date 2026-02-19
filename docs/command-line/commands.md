---
layout: default
title: BlobTools Commands
---

# BlobTools Commands

Reference guide for BlobTools command-line tools.

## Core Commands

### create

Create a new BlobDir dataset:

```bash
blobtools create --fasta assembly.fasta \
                 --meta metadata.json \
                 output_directory/
```

### add

Add data to an existing BlobDir:

```bash
blobtools add --hits blast_results.txt \
              --cov coverage.bam \
              output_directory/
```

### filter

Filter BlobDir sequences:

```bash
blobtools filter --query "gc > 0.4 AND read_cov > 10" \
                 --output_dir filtered/ \
                 input_directory/
```

### view

Export BlobDir data:

```bash
blobtools view --out tsv \
               input_directory/
```

## Additional Commands

- **validate** - Validate BlobDir structure
- **host** - Start interactive viewer
- **version** - Display version information

For detailed options, run:

```bash
blobtools [command] --help
```

## Next Steps

Check out [Tutorials](tutorials) for practical examples.
