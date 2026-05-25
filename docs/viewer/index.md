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

To visualize local BlobDir datasets, host the parent directory that contains one or more BlobDirs:

```bash
blobtools host --port 8080 /path/to/datasets/
```

Then open `http://localhost:8080/view/all` in your browser. See [Local Hosting](local-hosting) for the expected directory layout and Docker usage.
