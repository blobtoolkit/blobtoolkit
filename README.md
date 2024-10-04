# BlobToolKit (v4.3.13)

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/150091036.svg)](https://zenodo.org/badge/latestdoi/150091036)

BlobToolKit is described in our [BlobToolKit paper](https://doi.org/10.1534/g3.119.400908):

> BlobToolKit – Interactive quality assessment of genome assemblies
> Richard Challis, Edward Richards, Jeena Rajan, Guy Cochrane, Mark Blaxter
> G3: GENES, GENOMES, GENETICS April 1, 2020 vol. 10 no. 4 1361-1374;
> https://doi.org/10.1534/g3.119.400908

## About

Similar to [BlobTools v1](https://github.com/DRL/blobtools), **BlobTools2** is a command line tool designed to aid genome assembly QC and contaminant/cobiont detection and filtering. In addition to supporting interactive visualisation, a motivation for this reimplementation was to provide greater flexibility to include new types of information, such as [BUSCO](https://busco.ezlab.org) results and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) hit distributions.

**BlobToolKit** supports command-line filtering of datasets, assembly files and read files based on values or categories assigned to assembly contigs/scaffolds through the `blobtools filter` command. Interactive filters and selections made using the [BlobToolKit Viewer](https://github.com/blobtoolkit/viewer) can be reproduced on the command line and used to generate new, filtered datasets which retain all fields from the original dataset.

**BlobToolKit** is built around a file-based data structure, with data for each field contained in a separate `JSON` file within a directory (`BlobDir`) containing a single `meta.json` file with metadata for each field and the dataset as a whole. Additional fields can be added to an existing `BlobDir` using the `blobtools add` command, which parses an input to generate one or more additional `JSON` files and updates the dataset metadata. Fields are treated as generic datatypes, `Variable` (e.g. gc content, length and coverage), `Category` (e.g. taxonomic assignment based on BLAST hits) alongside `Array` and `MultiArray` datatypes to store information such as start, end, [NCBI](https://www.ncbi.nlm.nih.gov) taxid and bitscore for a set of blast hits to a single sequence. Support for new analyses can be added to **BlobTools2** by creating a new python module with an appropriate `parse` function.

To learn more about the development of the BlobTools approach, take a look at the papers by [Laetsch DR and Blaxter ML, 2017](https://f1000research.com/articles/6-1287/v1) and [Kumar et al., 2013](https://dx.doi.org/10.3389%2Ffgene.2013.00237).

## Installing

As of version 3.0.0, BlobTools2 and a local version of the BlobToolKit Viewer can be installed with:

```
pip install blobtoolkit
```

The `blobtools view`command requires firefox or a chromium-based browser to start the interactive viewer or to generate plots from the command line, these can be installed with:

```
conda install -c conda-forge firefox geckodriver
```

On MacOS, Xquartz is required to provide an X-windows environment. This can be installed without root sure privileges using:

```
brew install xquartz
```

## Example commands

Start the Viewer with an example dataset by using `_` as the dataset name (visit the URL shown in the command output)

```
blobtools view --local _
```

Generate a table of contigs in a filtered BlobDir:

```
blobtools filter --param length--Min=1000000 --table table.tsv /path/to/BlobDir
```

Generate a snail plot from a hosted BlobDir:

```
blobtools view --view snail --host https://blobtoolkit.genomehubs.org mSciVul1_1
```

To use a chromium-based browser (e.g. Google Chrome) in place of firefox, add `--driver chromium`

```
blobtools view --view snail --host https://blobtoolkit.genomehubs.org --driver chromium mSciVul1_1
```

See the [tutorial](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/) for details.

## Docker images

A set of Docker images are available from dockerhub:

- `genomehubs/blobtoolkit-blobtools` contains the blobtools executable, which includes the `blobtools view` command to run the Viewer and API
- `genomehubs/blobtoolkit` contains the blobtools executable along with all pipline dependencies
- `genomehubs/blobtoolkit-api` contains a standalone version of the BlobToolKit API
- `genomehubs/blobtoolkit-viewer` contains a standalone version of the BlobToolKit Viewer

The example commands above can all be run using the `genomehubs/blobtoolkit-blobtools` image. To access the API and Viewer running inside the container, it is necessary to map the default ports, 8000 and 8001:

```
docker run -it --rm --name blobtools -p 8000:8000 -p 8001:8001 genomehubs/blobtoolkit-blobtools:latest blobtools view --local _
```

## BlobToolKit pipeline

The BlobToolKit pipeline can be run by creating a YAML config file and environment variables to the `genomehubs/blobtoolkit` docker image.

```sh
# Set name of directory in which config.yaml file can be found
# The file should be available at /path/to/datasets/$ACCESSION/config.yaml
ACCESSION=GCA_963082805.1
# Set maximum number of threads to use
THREADS=16
# Set TRANSFER=true to remove intermediate files and place results in a separate directory
# at /path/to/output/<PREFIX>, where PREFIX is taken from assembly.prefix in config.yaml
# The final BlobDir will be available as /path/to/output/$PREFIX/PREFIX.tar
TRANSFER=true
docker run --rm \
    --name btk-$ACCESSION \
    -e ASSEMBLY=$ACCESSION \
    -e THREADS=$THREADS \
    -e TRANSFER=$TRANSFER \
    -v /path/to/datasets:/blobtoolkit/datasets \
    -v /path/to/databases:/blobtoolkit/databases \
    -v /path/to/output:/blobtoolkit/output \
    genomehubs/blobtoolkit:$RELEASE
```

Example config.yaml file:

```yaml
assembly:
  accession: GCA_963082805.1
  level: chromosome
  prefix: CAUJBB01
  scaffold-count: 197
  span: 377570513
busco:
  basal_lineages:
    - eukaryota_odb10
    - bacteria_odb10
    - archaea_odb10
  download_dir: /blobtoolkit/databases/busco_2021_06
  lineages:
    - endopterygota_odb10
    - insecta_odb10
    - arthropoda_odb10
    - metazoa_odb10
    - eukaryota_odb10
    - bacteria_odb10
    - archaea_odb10
reads:
  coverage:
    max: 30
  paired: []
  single:
    - base_count: 25744115650
      file: /blobtoolkit/datasets/GCA_963082805.1/reads/ERR11263500.fastq.gz
      platform: PACBIO_SMRT
      prefix: ERR11263500
settings:
  blast_chunk: 100000
  blast_max_chunks: 10
  blast_min_length: 1000
  blast_overlap: 0
  stats_chunk: 1000
  stats_windows:
    - 0.1
    - 0.01
    - 100000
    - 1000000
  taxdump: /blobtoolkit/databases/taxdump_2021_06
  tmp: /tmp
similarity:
  blastn:
    name: nt
    path: /blobtoolkit/databases/nt_2021_06
  defaults:
    evalue: 1.0e-10
    import_evalue: 1.0e-25
    max_target_seqs: 10
    taxrule: buscogenes
  diamond_blastp:
    import_max_target_seqs: 100000
    name: reference_proteomes
    path: /blobtoolkit/databases/uniprot_2021_06
    taxrule: blastp=buscogenes
  diamond_blastx:
    name: reference_proteomes
    path: /blobtoolkit/databases/uniprot_2021_06
taxon:
  name: Hemicrepidius niger
  taxid: "869179"
version: 1
```

## Contributing

If you find a problem or want to suggest a feature, please submit an issue.
