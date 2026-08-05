---
layout: default
title: BlobTools2
date: 2019-07-26
---

BlobTools2 is a reimplementation of [BlobTools](https://github.com/DRL/blobtools), written in [Python 3](https://www.python.org) with a fully modular design to make creating new datasets and adding additional analysis types even easier. Notable differences relative to BlobTools include the addition of a `blobtools add` command to allow new data to be added to an existing dataset and the ability to apply any filters from the interactive [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/) when using the`blobtools filter` and `blobtools view` commands.

```
BlobTools2 - assembly exploration, QC and filtering.

usage: blobtools [<command>] [<args>...] [-h|--help] [--version]

commands:
    add             add data to a BlobDir
    create          create a new BlobDir
    filter          filter a BlobDir
    host            host interactive view of all BlobDirs in a directory
    replace         call blobtools add with --replace flag
    view            generate plots using BlobToolKit Viewer
    -h, --help      show this
    -v, --version   show version number
See 'blobtools <command> --help' for more information on a specific command.
```

#### `blobtools create`

The minimum requirement to create a new dataset with BlobTools2 is an assembly FASTA file. The underlying data structure has been updated from a single JSON format BlobDB file to a collection of JSON files in a [BlobDir](https://blobtoolkit.genomehubs.org/specification/) directory. This change makes it trivial to separate dataset creation from the subsequent addition of analyses and presents the data in a format that can be efficiently processed for interactive visualisation with the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/).

#### `blobtools add`

Additional data can be added to an existing [BlobDir](https://blobtoolkit.genomehubs.org/specification/) by parsing analysis output files into one or more fields using the `blobtools add` command. This command can also be used to add metadata including links to external resources and full taxonomic information to a dataset. Currently supported analyses outputs include BLAST/Diamond sequence similarity searches, BAM/SAM/CRAM read mappings and BUSCO genome completeness assessments. Parsers are implemented as Python modules that convert the data to one of several generic datatypes (_identifier_, _variable_, _category_, _array_, _array of arrays_) so new analyses can be supported by adding an appropriate parser. The `blobtools replace` command calls `blobtools add` with a `--replace` flag to allow fields to be updated.

#### `blobtools filter`

Datasets can be filtered based on the values in any _variable_ or _category_ field, or using a list of _identifiers_. Filters may be applied to a complete dataset to allow for use of a reduced dataset without repeating analyses or applied to assembly FASTA and read FASTQ files to allow for reassembly and reanalysis. Filter parameters are all shared between BlobTools2 and the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/), allowing interactive sessions to be reproduced on the command line.

#### `blobtools host`

Assuming the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/) code and dependencies are available, the `blobtools host` command provides a convenient way to start the Viewer to begin interactive exploration of all [BlobDir](https://blobtoolkit.genomehubs.org/specification/) datasets in a directory.

#### `blobtools view`

Provides options to generate views by running the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/) from the command line.

See the _[BlobTools2 Tutorials](https://blobtoolkit.genomehubs.org/command-line/tutorials/)_ for more information on how to use BlobTools2 on your own datasets or check out our [open-source code](https://github.com/blobtoolkit/blobtools2) on GitHub to get started.
