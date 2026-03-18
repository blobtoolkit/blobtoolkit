---
layout: default
title: BlobDir Validator
---

# BlobDir Validator

The BlobDir validator ensures your datasets conform to the official BlobDir
specification. It is included with the `blobtools` command.

## Running the Validator

```bash
blobtools validate /path/to/blobdir/
```

You can also pass the `meta.json` file directly:

```bash
blobtools validate /path/to/blobdir/meta.json
```

## Validation Modes

By default the validator checks for all metadata fields required by the public
viewer instance, including assembly accession, prefix, taxon name and taxid.

To validate a **minimal BlobDir** (e.g. one created locally without a public
accession), pass `--basic` to disable those checks:

```bash
blobtools validate --basic /path/to/blobdir/
```

## Validation Output

### Success

```
VALID
```

Exit code 0.

### Failure

The validator prints an error message and exits with a non-zero code:

```
ERROR: /path/to/blobdir/meta.json
  - Missing required field: assembly.accession
  - Invalid value for gc: 1.23 (must be in range 0–1)
```

## What Is Checked

- Directory and file structure
- JSON schema compliance for `meta.json` and all field JSON files
- Data integrity (values within expected ranges)
- Required metadata fields (assembly accession, prefix, taxon name/id) unless
  `--basic` is used

## Dependencies

The validator requires `ujson` and `fastjsonschema`:

```bash
pip install ujson fastjsonschema
```

These are included when installing `blobtoolkit` via pip.

## Using Validated Datasets

A validated BlobDir can be safely used with:

- `blobtools host` — local viewer
- `blobtools filter` — command-line filtering and export
- BlobToolKit analysis pipelines
- Submission to the public BlobToolKit viewer

## Fixing Validation Errors

1. Read the error message carefully — the field name and issue are reported.
2. Refer to [BlobDir Format Specification](blobdir-format) for field
   definitions and allowed values.
3. Regenerate or repair the BlobDir (re-run `blobtools add` for the failing
   field, or manually correct `meta.json`).
4. Re-run `blobtools validate` to confirm the fix.
