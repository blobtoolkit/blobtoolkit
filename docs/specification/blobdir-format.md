---
layout: default
title: BlobDir Format Specification
---

# BlobDir Format Specification

Technical specification for the BlobDir JSON format used by BlobToolKit.

## Overview

BlobDir is a JSON-based format for storing genome assembly quality assessment data. It contains:

- Assembly sequence information
- Taxonomic assignments
- Coverage metrics
- Metadata

## Directory Structure

```
my_blobdir/
├── meta.json          # Metadata
├── blobdir.json       # Main data file
└── [optional files]
```

## meta.json Schema

```json
{
  "assembly_name": "string",
  "assembly_type": "string",
  "description": "string",
  "created": "ISO8601 timestamp",
  "version": "string"
}
```

## blobdir.json Schema

Main data structure containing:

- `records`: Array of contig records
- `categories`: Taxonomic categories
- `metrics`: Summary metrics

### Record Structure

```json
{
  "contig_id": "string",
  "length": "integer",
  "gc": "float",
  "coverage": "float",
  "taxonomy": "string",
  "hits": []
}
```

## JSON Schema Validation

The complete JSON schema is available in the [Validator](validator.md).

## Compatibility

Current version: **1.0**

Older versions may have different schemas. Use the validator to check compatibility.
