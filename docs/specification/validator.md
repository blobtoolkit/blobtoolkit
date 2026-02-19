---
layout: default
title: BlobDir Validator
---

# BlobDir Validator

The BlobDir validator ensures your datasets conform to the official BlobDir specification.

## Running the Validator

```bash
blobtools validate /path/to/blobdir/
```

## Validation Output

The validator checks:

- Directory structure
- JSON schema compliance
- Data integrity
- Required fields
- Data type validation

### Success Output

```
✓ BlobDir is valid
Version: 1.0
Records: 1234
```

### Error Output

The validator provides detailed error messages:

```
✗ Validation failed:
  - Missing required field: meta.json
  - Invalid JSON in blobdir.json: line 42
  - Record 'contig_123': invalid gc value (should be 0-1)
```

## Using Validated Data

Validated BlobDirs can be safely used with:

- BlobTools viewers
- Analysis pipelines
- Export functions

## Fixing Validation Errors

1. Check error messages carefully
2. Refer to [BlobDir Format Specification](blobdir-format)
3. Regenerate or repair the BlobDir
4. Re-run validator to confirm

## Custom Validation

To check specific aspects:

```bash
blobtools validate --check schema /path/to/blobdir/
blobtools validate --check structure /path/to/blobdir/
```
