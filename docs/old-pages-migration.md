---
layout: default
title: Old Pages Migration Audit
---

# Old Pages Migration Audit

This page tracks migration from `docs/old_pages` into the current docs structure.

## Command Truth Used for This Audit

Current CLI support is derived from source entry points:

- `blobtools` subcommands: `add`, `create`, `filter`, `host`, `remove`, `replace`, `validate`, `view`
- `btk` subcommands: `pipeline`, `blobtools`

Preferred pipeline user guidance should point to:

- https://pipelines.tol.sanger.ac.uk/blobtoolkit

Current status: preferred Nextflow pipeline is configured and run outside the BlobToolKit CLI.

These were validated against:

- `setup.py` entry points
- `src/blobtools/blobtools.py`
- `src/btk/btk.py`

## Folder-Level Summary

- `old_pages/blobtools2`: mostly still relevant, migrate into `docs/command-line`
- `old_pages/btk-viewer`: mostly still relevant, migrate into `docs/viewer`
- `old_pages/pipeline`: mostly relevant but contains legacy Snakemake details, migrate into `docs/pipeline` with caveats
- `old_pages/specification`: relevant, merge into `docs/specification`
- `old_pages/install`: relevant but outdated versions, migrate carefully into `docs/installation`
- `old_pages/home/about/tools`: mostly historical/context pages, optional migration

## Migration Matrix

- `old_pages/blobtools2/blobtools2-tutorials/creating-a-dataset/index.md` -> `docs/command-line/tutorials.md` -> MERGE
- `old_pages/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/index.md` -> `docs/command-line/tutorials.md` -> MERGE
- `old_pages/blobtools2/blobtools2-tutorials/filtering-a-dataset/index.md` -> `docs/command-line/tutorials.md` and `docs/command-line/commands.md` -> MERGE
- `old_pages/pipeline/pipeline-tutorials/configuring-the-pipeline/index.md` -> `docs/pipeline/configuration.md` -> MERGE
- `old_pages/install/index.md` (database section) -> `docs/installation/databases.md` -> SPLIT + MERGE
- `old_pages/specification/specification-tutorials/validating-datasets/index.md` -> `docs/specification/validator.md` -> DONE
- `old_pages/btk-viewer/viewer-tutorials/hosting-a-local-instance/index.md` -> `docs/viewer/local-hosting.md` -> DONE
- `old_pages/btk-viewer/viewer-tutorials/exploring-views/index.md` -> `docs/viewer/visualisations.md` -> DONE
- `old_pages/btk-viewer/viewer-tutorials/adjusting-plot-settings/index.md` -> `docs/viewer/visualisations.md` -> DONE (merged)
- `old_pages/btk-viewer/viewer-tutorials/filtering-assemblies/index.md` -> `docs/viewer/visualisations.md` -> DONE (merged, ## Filtering Assemblies)
- `old_pages/btk-viewer/viewer-tutorials/changing-plot-axes/index.md` -> `docs/viewer/visualisations.md` -> DONE (merged, ## Changing Plot Axes)
- `old_pages/btk-viewer/viewer-tutorials/using-selections/index.md` -> `docs/viewer/visualisations.md` -> DONE (merged, ## Using Selections)
- `old_pages/btk-viewer/viewer-tutorials/viewing-large-datasets/index.md` -> `docs/viewer/visualisations.md` -> DONE (merged, ## Performance subsections)
- `old_pages/btk-viewer/viewer-tutorials/public-vs-local-instances/index.md` -> `docs/viewer/local-hosting.md` -> DONE (merged, ## Public vs. Local Instances)
- `old_pages/btk-viewer/viewer-tutorials/searching-available-datasets/index.md` -> `docs/viewer/local-hosting.md` -> DONE (merged, ## Searching Available Datasets)
- `old_pages/pipeline/pipeline-tutorials/running-the-pipeline-on-a-cluster/index.md` -> `docs/pipeline/advanced/cluster-deployment.md` -> DONE
- `old_pages/pipeline/pipeline-tutorials/running-the-pipeline-in-a-container/index.md` -> `docs/pipeline/advanced/container-deployment.md` -> DONE

## Batch Plan

### Batch 1 (Completed)

- Expanded command-line docs to include migrated workflows and current command set.
- Expanded pipeline getting-started and configuration docs.
- Added installation database setup page.
- Updated navigation for new installation page.

### Batch 2 (Completed)

- Expanded viewer docs: local-hosting with environment variable table, port
  forwarding details and performance notes.
- Expanded viewer visualisations with all view types, settings reference table
  and performance guidance (merged adjusting-plot-settings content).
- Expanded specification validator with --basic mode, dependencies, and fix
  workflow.

### Batch 3 (Completed)

- Added `docs/pipeline/advanced/cluster-deployment.md`: Nextflow (primary) +
  legacy Snakemake cluster job submission, drmaa per-rule job example.
- Added `docs/pipeline/advanced/container-deployment.md`: Nextflow profiles
  (primary) + legacy Docker manual run with volume mount guide.
- Updated navigation to include both new advanced pages.

### Batch 4 (Completed)

- Extended `docs/viewer/visualisations.md` with:
  - ## Performance: Plot graphics threshold, Circle limit, Static threshold subsections with images.
  - ## Changing Plot Axes: cov-cov plots, category rank switching with images.
  - ## Filtering Assemblies: slider, category-toggle, selection-based, invert with images.
  - ## Using Selections: BUSCO-based selection, tracking across views, table inspection with images.
- Extended `docs/viewer/local-hosting.md` with:
  - ## Public vs. Local Instances: home page differences, search results table with images.
  - ## Searching Available Datasets: search box autocomplete, sort/customise, taxonomy browser with images.
- Copied all six remaining viewer tutorial image sets to `docs/assets/img/viewer/` (138 images total).

## Contradiction Checks to Keep Applying

- Replace old `BlobTools2` naming with `blobtools`.
- Prefer https://pipelines.tol.sanger.ac.uk/blobtoolkit for pipeline config/run guidance.
- Keep `btk pipeline` references only where explicitly labeled as legacy/local.
- Remove outdated Python and Node version references.
- Avoid linking to deleted `old_pages` URLs in migrated content.
