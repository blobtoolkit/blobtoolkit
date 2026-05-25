---
layout: default
title: Local Hosting
---

# Local Hosting with BlobTools Host

The BlobTools host command starts an interactive web server for visualizing your blob plot data locally. The command expects a directory containing one or more BlobDirs.

## Starting the Viewer

For a single BlobDir, pass the parent directory that contains the BlobDir:

```bash
blobtools host --port 8080 /path/to/datasets/
```

where the directory contains a BlobDir such as:

```text
datasets/
`-- example-1/
    `-- meta.json
```

Then open your browser and navigate to `http://localhost:8080/view/all`. You can also search for `all` in the viewer to list available datasets.

## Hosting Multiple BlobDirs

To host more than one local dataset, place each BlobDir in the same parent directory and host that parent:

```text
datasets/
|-- example-1/
|   `-- meta.json
`-- example-2/
    `-- meta.json
```

```bash
blobtools host --api-port 8000 --port 8080 /path/to/datasets/
```

Open `http://localhost:8080/view/all` to browse the indexed datasets, or open a dataset directly using a URL such as `http://localhost:8080/view/all/dataset/example-1`.

With Docker, mount the parent directory to `/blobtoolkit/datasets` and host that mounted path:

```bash
docker run -d --rm --name btk \
    -v "$PWD/datasets:/blobtoolkit/datasets" \
    -p 8000:8000 -p 8080:8080 \
    genomehubs/blobtoolkit:latest \
    blobtools host --api-port 8000 --port 8080 /blobtoolkit/datasets
```

## Port Configuration

To use a different port:

```bash
blobtools host --port 3000 /path/to/datasets/
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
