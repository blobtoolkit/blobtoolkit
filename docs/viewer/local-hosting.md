---
layout: default
title: Local Hosting
---

# Local Hosting with BlobTools Host

The BlobTools host command starts an interactive web server for visualizing your blob plot data locally.

## Starting the Viewer

```bash
blobtools host --port 8080 /path/to/blobdir/
```

Then open your browser and navigate to `http://localhost:8080`

## Port Configuration

To use a different port:

```bash
blobtools host --port 3000 /path/to/blobdir/
```

## Connecting from Remote Servers

If BlobTools is running on a remote server:

```bash
ssh -L 8080:localhost:8080 user@remote-server
# Then access http://localhost:8080 locally
```

## Features

The web viewer provides:

- Interactive blob plots
- Contig filtering and selection
- Taxonomic distribution visualization
- Data export capabilities

## Next Steps

Learn how to interpret the visualizations in [Visualisations Guide](visualisations).
