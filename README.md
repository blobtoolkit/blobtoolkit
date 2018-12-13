# BlobTools 2.0

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Build Status](https://travis-ci.org/blobtoolkit/blobtools2.svg?branch=master)](https://travis-ci.org/blobtoolkit/blobtools-add)
[![Coverage Status](https://coveralls.io/repos/github/blobtoolkit/blobtools-add/badge.svg?branch=master)](https://coveralls.io/github/blobtoolkit/blobtools2?branch=master)

A new implementation of [BlobTools](https://github.com/DRL/blobtools) with support for interactive data exploration via the [BlobToolKit Viewer](https://github.com/blobtoolkit/viewer).

## About

Similar to [BlobTools v1](https://github.com/DRL/blobtools), **BlobTools2** is a command line tool designed to aid genome assembly QC and contaminant/cobiont detection and filtering. In addition to supporting interactive visualisation, a motivation for this reimplementation was to provide greater flexibility to include new types of information, such as [BUSCO](https://busco.ezlab.org) results and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) hit distributions.

**BlobTools2** supports command-line filtering of datasets, assembly files and read files based on values or categories assigned to assembly contigs/scaffolds through the `blobtools filter` command. Interactive filters and selections made using the [BlobToolKit Viewer](https://github.com/blobtoolkit/viewer) can be reproduced on the command line and used to generate new, filtered datasets which retain all fields from the original dataset.

**BlobTools2** is built around a file-based data structure, with data for each field contained in a separate `JSON` file within a directory (`BlobDir`) containing a single `meta.json` file with metadata for each field and the dataset as a whole. Additional fields can be added to an existing `BlobDir` using the `blobtools add` command, which parses an input to generate one or more additional `JSON` files and updates the dataset metadata. Fields are treated as generic datatypes, `Variable` (e.g. gc content, length and coverage), `Category` (e.g. taxonomic assignment based on BLAST hits) alongside `Array` and `MultiArray` datatypes to store information such as start, end, [NCBI](https://www.ncbi.nlm.nih.gov) taxid and bitscore for a set of blast hits to a single sequence. Support for new analyses can added to **BlobTools2** by creating a new python module with an appropriate `parse` function


## Installing

These install instructions assume a Unix/Linux system with standard development tools and [Conda](https://conda.io/docs/user-guide/install/index.html) installed. The [BlobToolKit Viewer](https://github.com/blobtoolkit/viewer) (which is included as a submodule) also requires `libpng-dev`.

1. Create and activate a Conda environment
```
conda create -n blobtools2 -y python=3.6 docopt pyyaml ujson pysam tqdm nodejs seqtk
conda activate blobtools2
```

2. Download and extract the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) taxdump:
```
mkdir -p /path/to/new_taxdump
cd /path/to/new_taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
cd -
```

3. Clone this repository:
```
git clone --recursive https://github.com/blobtoolkit/blobtools2
```

4. Install node modules for [BlobToolKit Viewer](https://github.com/blobtoolkit/viewer):
```
cd blobtools2/viewer
npm install
cd -
```


## Examples

```
./blobtools create --fasta examples/assembly.fasta tmp/example
```

```
./blobtools add --cov examples/assembly.reads.bam tmp/example
```

```
./blobtools add --hits examples/blast.out --hits examples/diamond.out --taxdump /path/to/new_taxdump/ --taxrule bestsumorder tmp/example
```

```
./blobtools add --busco examples/busco.tsv tmp/example
```


```
./blobtools create --meta meta.yaml --fasta assembly.fasta --key record_type=scaffold path/to/dataset
```

```
./blobtools add --cov mapping_file.bam --taxdump path/to/new_taxdump --hits blast.out --threads 8 path/to/dataset
```

```
./blobtools filter --param gc--Max=0.26 --output path/to/filtered_dataset path/to/dataset
```

```
./blobtools add --key plot.x=gc --key plot.y=mapping_file_cov --key plot.z=length --key plot.cat=bestsum_phylum path/to/filtered_dataset
```


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
