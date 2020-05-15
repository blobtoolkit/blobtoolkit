# BlobTools2

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Build Status](https://travis-ci.org/blobtoolkit/blobtools2.svg?branch=master)](https://travis-ci.org/blobtoolkit/blobtools2)
[![Coverage Status](https://coveralls.io/repos/github/blobtoolkit/blobtools2/badge.svg?branch=master)](https://coveralls.io/github/blobtoolkit/blobtools2?branch=master)
[![DOI](https://zenodo.org/badge/150091036.svg)](https://zenodo.org/badge/latestdoi/150091036)

A new implementation of [BlobTools](https://github.com/DRL/blobtools) with support for interactive data exploration via the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer).

More information and tutorials are available at [blobtoolkit.genomehubs.org/blobtools2/](https://blobtoolkit.genomehubs.org/blobtools2/)

BlobToolKit is described in our [BlobToolKit paper](https://doi.org/10.1534/g3.119.400908):

> BlobToolKit – Interactive quality assessment of genome assemblies
> Richard Challis, Edward Richards, Jeena Rajan, Guy Cochrane, Mark Blaxter
> G3: GENES, GENOMES, GENETICS April 1, 2020 vol. 10 no. 4 1361-1374;
> https://doi.org/10.1534/g3.119.400908

## About

Similar to [BlobTools v1](https://github.com/DRL/blobtools), **BlobTools2** is a command line tool designed to aid genome assembly QC and contaminant/cobiont detection and filtering. In addition to supporting interactive visualisation, a motivation for this reimplementation was to provide greater flexibility to include new types of information, such as [BUSCO](https://busco.ezlab.org) results and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) hit distributions.

**BlobTools2** supports command-line filtering of datasets, assembly files and read files based on values or categories assigned to assembly contigs/scaffolds through the `blobtools filter` command. Interactive filters and selections made using the [BlobToolKit Viewer](https://github.com/blobtoolkit/viewer) can be reproduced on the command line and used to generate new, filtered datasets which retain all fields from the original dataset.

**BlobTools2** is built around a file-based data structure, with data for each field contained in a separate `JSON` file within a directory (`BlobDir`) containing a single `meta.json` file with metadata for each field and the dataset as a whole. Additional fields can be added to an existing `BlobDir` using the `blobtools add` command, which parses an input to generate one or more additional `JSON` files and updates the dataset metadata. Fields are treated as generic datatypes, `Variable` (e.g. gc content, length and coverage), `Category` (e.g. taxonomic assignment based on BLAST hits) alongside `Array` and `MultiArray` datatypes to store information such as start, end, [NCBI](https://www.ncbi.nlm.nih.gov) taxid and bitscore for a set of blast hits to a single sequence. Support for new analyses can be added to **BlobTools2** by creating a new python module with an appropriate `parse` function.

To learn more about the development of the BlobTools approach, take a look at the papers by [Laetsch DR and Blaxter ML, 2017](https://f1000research.com/articles/6-1287/v1) and [Kumar et al., 2013](https://dx.doi.org/10.3389%2Ffgene.2013.00237).


## Installing

The most up to date instructions can be found at [blobtoolkit.genomehubs.org/install/](https://blobtoolkit.genomehubs.org/install/).

## Examples

### Create a new dataset
- [blobtools create](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/creating-a-dataset/)

### Adding more data

- [Coverage](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/adding-coverage/)
- [BLAST/Diamond hits](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/adding-hits/)
- [BUSCO results](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/adding-busco/)

### Setting dataset metadata

- [File-based metadata](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/creating-a-dataset/)
- [Adding taxonomy information](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/updating-metadata/)
- [Updating individual keys](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/updating-metadata/)
- [Adding external links](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/updating-metadata/)

### Filtering datasets

- [Specifying filter parameters](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/filtering-a-dataset/)
- [Filtering data files](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/filtering-a-dataset/)

### Visualising datasets

- [Interactive visualisation](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/opening-a-dataset-in-the-viewer/)
- [Command line visualisation](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/generating-plots-on-the-command-line/)


## Contributing

If you find a problem or want to suggest a feature, please submit an issue.

If you want to contribute code, pull requests are welcome but please make sure your code passes the linting, testing and style checks before submitting a pull request. Additional development dependencies are listed in `requirements.txt`

Run linter/testing locally:
```
./run_tests.sh
```

Set up pre-commit hook to automate test running:
```
ln -s ../../pre-commit.sh .git/hooks/pre-commit
```
