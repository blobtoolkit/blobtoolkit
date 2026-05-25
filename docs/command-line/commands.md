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

Export BlobDir data or plots:

```bash
blobtools view --out tsv \
               input_directory/
```

Blob plots use the category field stored in the BlobDir plot metadata by
default. Hits added with the default taxrule create category fields named
`bestsumorder_<rank>`, such as `bestsumorder_phylum`,
`bestsumorder_family`, and `bestsumorder_genus`.

Use the `catField` parameter to render a blob plot at a different taxonomic
rank from the command line:

```bash
blobtools view --plot \
               --view blob \
               --format png \
               --out ./ \
               --param catField=bestsumorder_family \
               input_directory/
```

The same parameter can be passed without `--plot` when rendering through the
viewer-backed export path:

```bash
blobtools view --view blob \
               --format png \
               --out ./ \
               --param catField=bestsumorder_genus \
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
