---
layout: default
title: Viewer
nav_order: 5
---

# BlobToolKit Viewer

Interactive web-based visualization of genome assembly quality and taxonomic composition.

## Overview

This section covers:

- **[Local Hosting](local-hosting)** - Running the viewer on your machine
- **[Visualisations](visualisations)** - Interpreting plots and metrics

## Key Features

The viewer provides:

- Interactive blob plots with real-time filtering
- Snail plots showing cumulative contig coverage
- Taxonomic distribution analysis
- Contig-level data exploration
- Export capabilities for publication-quality plots

## Quick Start

To visualize a local BlobDir dataset:

```bash
blobtools host --port 8080 /path/to/blobdir/
```

Then open `http://localhost:8080` in your browser.
