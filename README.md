# BlobTools2

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Build Status](https://travis-ci.org/blobtoolkit/blobtools2.svg?branch=master)](https://travis-ci.org/blobtoolkit/blobtools2)
[![Coverage Status](https://coveralls.io/repos/github/blobtoolkit/blobtools2/badge.svg?branch=master)](https://coveralls.io/github/blobtoolkit/blobtools2?branch=master)
[![DOI](https://zenodo.org/badge/150091036.svg)](https://zenodo.org/badge/latestdoi/150091036)


A new implementation of [BlobTools](https://github.com/DRL/blobtools) with support for interactive data exploration via the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer).


## About

Similar to [BlobTools v1](https://github.com/DRL/blobtools), **BlobTools2** is a command line tool designed to aid genome assembly QC and contaminant/cobiont detection and filtering. In addition to supporting interactive visualisation, a motivation for this reimplementation was to provide greater flexibility to include new types of information, such as [BUSCO](https://busco.ezlab.org) results and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) hit distributions.

**BlobTools2** supports command-line filtering of datasets, assembly files and read files based on values or categories assigned to assembly contigs/scaffolds through the `blobtools filter` command. Interactive filters and selections made using the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer) can be reproduced on the command line and used to generate new, filtered datasets which retain all fields from the original dataset.

**BlobTools2** is built around a file-based data structure, with data for each field contained in a separate `JSON` file within a directory (`BlobDir`) containing a single `meta.json` file with metadata for each field and the dataset as a whole. Additional fields can be added to an existing `BlobDir` using the `blobtools add` command, which parses an input to generate one or more additional `JSON` files and updates the dataset metadata. Fields are treated as generic datatypes, `Variable` (e.g. gc content, length and coverage), `Category` (e.g. taxonomic assignment based on BLAST hits) alongside `Array` and `MultiArray` datatypes to store information such as start, end, [NCBI](https://www.ncbi.nlm.nih.gov) taxid and bitscore for a set of blast hits to a single sequence. Support for new analyses can be added to **BlobTools2** by creating a new python module with an appropriate `parse` function.

To learn more about the BlobTools approach, take a look at the papers by [Laetsch DR and Blaxter ML, 2017](https://f1000research.com/articles/6-1287/v1) and [Kumar et al., 2013](https://dx.doi.org/10.3389%2Ffgene.2013.00237). **BlobTools2** is intended to to become a more flexible replacement for [BlobTools v1](https://github.com/DRL/blobtools) but not all of the functionality has been added yet. If you need functions that are not yet ready in **BlobTools2**, please continue to use (and cite!) Dom Laetsch's [previous version](https://github.com/DRL/blobtools) and submit a feature request on the [issue tracker](https://github.com/blobtoolkit/blobtools2/issues) if you'd like to see it added.


## Installing

These install instructions assume a Unix/Linux system with standard development tools and [Conda](https://conda.io/docs/user-guide/install/index.html) installed. Installing the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer) as described below also requires `libpng-dev`.

1. Create and activate a Conda environment
```
conda create -n blobtools2 -y python=3.6 docopt pyyaml ujson tqdm nodejs seqtk
conda activate blobtools2
conda install -c conda-forge geckodriver selenium pysam
```

2. Download and extract the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) taxdump:
```
mkdir -p taxdump
cd taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
cd ..
```

3. Clone this repository:
```
git clone https://github.com/blobtoolkit/blobtools2
```

4. Clone and install the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer):
```
git clone https://github.com/blobtoolkit/viewer
cd viewer
npm install
cd ..
```


## Examples

### Create a new dataset

Either add all data in a single command to generate a full dataset ready to explore in the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer):
```
./blobtools create --fasta examples/assembly.fasta --cov examples/assembly.reads.bam --hits examples/blast.out --taxdump ../taxdump tmp/dataset_1
```

Or create a new dataset based on an assembly and add further fields later:
```
./blobtools create --fasta examples/assembly.fasta tmp/dataset_2
```

### Import an existing dataset from a BlobDB file

If you already have a `blobDB.json` file from [BlobTools v1](https://github.com/DRL/blobtools), this can be converted to the BlobTools2 `BlobDir` format:
```
./blobtools create --blobdb examples/blobDB.json tmp/dataset_3
```

### Adding more data

All `blobtools add` flags can also be used with `blobtools create` to generate a more complete dataset in a single command.

#### Coverage

Add coverage information from BAM or CRAM files:
```
./blobtools add --cov examples/assembly.reads.bam tmp/dataset_2
```

Specify `--pileup-args` to customise coverage calculations (e.g. as the example dataset is so small setting `stepper=nofilter` ensures all aligned positions are counted):
```
./blobtools add --cov examples/assembly.reads.bam --pileup-args stepper=nofilter tmp/dataset_2
```

To speed up coverage file processing, run commands using multiple threads:
```
./blobtools add --cov examples/assembly.reads.bam  --threads 16 tmp/dataset_2
```

Field names are based on coverage filenames, to explicitly set a field name, add `=<fieldname>` immediately after the coverage filename:
```
./blobtools add --cov examples/assembly.reads.bam=library1 --threads 16 tmp/dataset_2
```

#### BLAST/Diamond hits

Sequence similarity search results are used to assign taxonomy based on [BlobTools v1](https://github.com/DRL/blobtools) taxrules. In order to process taxonomic information, a local copy of the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) taxdump must be available:

```
./blobtools add --hits examples/blast.out --hits examples/diamond.out --taxdump ../taxdump --taxrule bestsumorder tmp/dataset_2
```

#### BUSCO results

Results from comparison against one or more [BUSCO](https://busco.ezlab.org) sets can be imported:
```
./blobtools add --busco examples/busco.tsv tmp/dataset_2
```

### Setting dataset metadata

#### File-based metadata
Metadata can be loaded from a `YAML` or `JSON` format file. Any fields in an `assembly` or `taxon` section will be indexed by the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer) API and will be searchable in the viewer:
```
record_type: scaffold
assembly:
  accession: GCA_000950515.2
  alias: O_ochengi_Ngaoundere
  bioproject: PRJEB1204
  biosample: SAMEA1034766
  prefix: FJNM01
taxon:
  name: Onchocerca ochengi
  taxid: 42157
  phylum: Nematoda
```

```
./blobtools add --meta examples/meta.yaml tmp/dataset_2
```

#### Taxonomic ranks
Full taxonomic lineage can be loaded for a given taxid based on information in the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) taxdump:
```
./blobtools add --taxid 42157 --taxdump ../blobtools-add/taxdump tmp/dataset_2
```

#### Individual keys
Specific keys in the metadata can be edited directly using the `--key` flag:
```
./blobtools add --key taxon.taxid=42157 --key taxon.name="Onchocerca ochengi" --key ./assembly.accession=draft tmp/dataset_2
```

#### External links
Links to external resources can be added using the `--link` flag, these will be shown alongside the appropriate data in the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer). The location within the metadata is specified by a `.` delimited string with the last part being a title for the link. The url is parsed to replace `{key}` with the corresponding value from the same metadata subsection:
```
./blobtools add --link taxon.taxid.ENA="https://www.ebi.ac.uk/ena/data/view/Taxon:{taxid}" --link taxon.name.Wikipedia="https://en.wikipedia.org/wiki/{name}" tmp/dataset_2
```
By default **BlobTools2** will test that the link target exists before adding the link to the metadata. to disable this behaviour, use the `--skip-link-test` flag:
```
./blobtools add --link taxon.name.Wikipedia="https://en.wikipedia.org/wiki/{name}" --skip-link-test tmp/dataset_2
```

Links can also be added to individual records:
```
--link record.ENA="https://www.ebi.ac.uk/ena/data/view/{id}" tmp/dataset_2
```

and to BLAST hits:
```
./blobtools add --link position.NCBI="https://www.ncbi.nlm.nih.gov/nuccore/{subject}" tmp/dataset_2
```

To link to different resources for hits from different files, links can be specified with an index used in the same order as the files were listed for `blobtools add --hits`:
```
./blobtools add --link position.0.NCBI="https://www.ncbi.nlm.nih.gov/nuccore/{subject}" --link position.1.UniProt="https://www.uniprot.org/uniprot/{subject}" tmp/dataset_2
```

### Filtering datasets

Datasets can be filtered interactively using the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer) or directly on the command line. Most options in the viewer are captured in the URL so interactive filtering can be reproduced on the command line from the values in the URL query string.

#### Specifying filter parameters
Variable based filters can be specified individually. Use the `--output` flag to specify an output directory:
```
./blobtools filter --param length--Min=5000 --output tmp/example_len_gt_5000 tmp/dataset_2
```

Or by pasting a query string or complete URL:
```
./blobtools filter --query-string "http://localhost:8080/all/dataset/example/blob?length--Min=5000#Filters" --output tmp/example_len_gt_5000 tmp/dataset_2
```

Available filters for `Variable` fields are:
```
<field_id>--Min - lowest value to include
<field_id>--Max - highest value to include
<field_id>--Inv - include values outside the range specified by --Min and --Max
```

Category filters operate on keys:
```
./blobtools filter --param bestsum_phylum--Keys=no-hit --output tmp/example_len_gt_5000 tmp/dataset_2
```

Available filters for `Category` fields are:
```
<field_id>--Keys - comma-separated list of strings matching category names or integers matching category keys to exclude
<field_id>--Inv - include rather than exclude --Keys
```

Selection-based filters are not captured in the query string so to reproduce an interactive selection on the command line, it is necessary to export the current selection from the viewer as a list, which can be loaded using the `--json` flag:
```
./blobtools filter --json examples/list.json --output tmp/example_len_gt_5000 tmp/dataset_2
```

All filters can be inverted to make them inclusive rather than exclusive:
```
./blobtools filter --json examples/list.json --invert --output tmp/example_len_lt_5000 tmp/dataset_2
```

#### Filtering data files
As an alternative to (or in addition to) creating a new, filtered dataset, filters can also be used to obtain subsets of assembly and read files based on the identifiers in the filtered set:

```
./blobtools filter --param length--Min=5000 --fasta examples/assembly.fasta tmp/dataset_2
```

```
./blobtools filter --param gc--Max=0.26 --fastq examples/reads_1.fq.gz --fastq examples/reads_2.fq.gz --cov examples/assembly.reads.bam  tmp/dataset_2
```

### Visualising datasets

#### Interactive visualisation
Datasets can be visualised interactively by running the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer). the `host` command provides a convenient way to run a local instance of the viewer hosting all datasets in the specified directory:
```
./blobtools host tmp
```

Note that this runs code in a non-optimised development mode not suited to public hosting. To make an instance of the viewer available publicly, see the more detailed instruction in the viewer [repository](https://github.com/blobtoolkit/viewer).

By default this starts the API on port 8000 and the viewer on port 8080, to change these settings pass additional flags as below:
```
./blobtools host --port 8081 --api-port 8001 tmp
```

Optionally specify the `--hostname` flag to allow connections from other hosts:
```
./blobtools host --port 8081 --api-port 8001 --hostname $(hostname) tmp
```

#### Command line visualisation
This functionality will be added to **BlobTools2** but for now, the only way to generate images on the command line is to use the `cli` utility within the [BlobToolKit viewer](https://github.com/blobtoolkit/viewer). This command line interface runs a [Firefox](https://www.mozilla.org/en-GB/firefox/new/) browser in headless mode to allow `png` and `svg` images of the main plot types to be directly generated with any query string options applied.


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
