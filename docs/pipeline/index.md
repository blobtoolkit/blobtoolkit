---
layout: default
title: Pipeline
nav_order: 3
---

# BlobToolKit Pipeline

The BlobToolKit pipeline automates genome assembly quality assessment end-to-end, from data retrieval through analysis and visualization.

## Preferred Pipeline

The preferred pipeline for configuration and execution is the Nextflow implementation documented at:

- https://pipelines.tol.sanger.ac.uk/blobtoolkit

At present this preferred Nextflow workflow is run outside the BlobToolKit CLI.

## Quick Links

<!-- QUICK_LINKS -->

## Overview

This section covers:

- **[Getting Started](getting-started)** - preferred workflow and legacy local options
- **[Configuration](configuration)** - legacy/local configuration reference

## Pipeline Components

The pipeline includes:

- Automated read mapping and coverage calculation
- BUSCO gene annotation and analysis
- Taxonomic assignment via BLAST/Diamond
- Data aggregation and BlobDir generation

## Running the Pipeline

For the recommended approach, use the external pipeline documentation at https://pipelines.tol.sanger.ac.uk/blobtoolkit.

Use the pages in this section for context and legacy/local pipeline notes.
