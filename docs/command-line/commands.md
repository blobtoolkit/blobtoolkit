---
layout: default
title: BlobTools Commands
---

# BlobTools Commands

Reference guide for BlobTools command-line tools.

The supported subcommands are defined in the Python entry points and include:

- `add`
- `create`
- `filter`
- `host`
- `remove`
- `replace`
- `validate`
- `view`

## Core Commands

### create

Create a new BlobDir dataset. `create` is an alias of `add` with creation mode:

```bash
blobtools create --fasta assembly.fasta \
                 --meta metadata.json \
                 output_directory/
```

### add

Add one or more analyses to an existing BlobDir:

```bash
blobtools add --hits blast_results.txt \
              --cov coverage.bam \
              output_directory/
```

### filter

Filter BlobDir sequences using numeric ranges or category keys:

```bash
blobtools filter --query "gc > 0.4 AND read_cov > 10" \
                 --output_dir filtered/ \
                 input_directory/
```

Common filter options:

- `--param` to set one filter parameter at a time, e.g. `length--Min=1000`
- `--query-string` to reuse parameters copied from a viewer URL
- `--json` to load a selection/list exported from the viewer
- `--invert` to invert matching behavior
- `--fasta`, `--fastq`, `--text` to filter external files
- `--summary`, `--table` for summary or tabular outputs

### view

Export BlobDir data:

```bash
blobtools view --out tsv \
               input_directory/
```

### validate

Validate BlobDir structure and content:

```bash
blobtools validate /path/to/blobdir/
```

### host

Run a local server for interactive exploration:

```bash
blobtools host --port 8080 /path/to/blobdir/
```

### remove

Remove one or more fields from a BlobDir:

```bash
blobtools remove --fields bestsumorder_phylum /path/to/blobdir/
```

### replace

Replace existing values when re-importing data. Internally this maps to `add --replace`:

```bash
blobtools replace --hits updated_hits.tsv /path/to/blobdir/
```

## Additional Commands

- **version** - Display version information (`blobtools --version`)

For detailed options, run:

```bash
blobtools [command] --help
```

## Next Steps

Check out [Tutorials](tutorials) for practical examples.
