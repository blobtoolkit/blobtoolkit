---
title: Creating a dataset
date: 2020-02-10
---

The minimum requirement to create a new dataset with BlobTools2 is an assembly FASTA file, but it is possible to create a much richer dataset by adding additional data and metadata. The underlying [BlobDir](https://blobtoolkit.genomehubs.org/specification/) data structure makes it trivial to separate dataset creation from the subsequent addition of analyses (see _[Adding data to a dataset](https://blobtoolkit.genomehubs.org/command-line/tutorials/adding-data-to-a-dataset/)_) and presents the data in a format that can be efficiently processed for interactive visualisation with the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/).

To follow along with this tutorial using your own files, just change the names of any files below written in `ALL_CAPS`.

#### Prepare files

Files for this example should all be saved in the same directory (`~/BTK_TUTORIAL/FILES`):

mkdir -p ~/BTK_TUTORIAL/FILES
cd ~/BTK_TUTORIAL/FILES

Make sure you have a genome sequence file (`ASSEMBLY_NAME.fasta`). These examples use a publicly available _Strongyloides venezuelensis_ assembly:

```
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/028/725/GCA_001028725.1_S_venezuelensis_HH1/GCA_001028725.1_S_venezuelensis_HH1_genomic.fna.gz | gunzip -c > ASSEMBLY_NAME.fasta
```

(Optional) Create a metadata file (`ASSEMBLY_NAME.yaml`) in the same directory as the assembly. As well as providing a way to link metadata to a dataset, entries in the assembly and taxon sections of the metadata are searchable within the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/). The entries in this file overlap with the full metadata file used by the [BlobToolKit Pipeline](https://blobtoolkit.genomehubs.org/pipeline/) (see [_configuring the pipeline_](https://blobtoolkit.genomehubs.org/pipeline/pipeline-tutorials/configuring-the-pipeline/)):

```
nano ASSEMBLY_NAME.yaml
```

```
assembly:
  accession: GCA_001028725.1
  alias: S_venezuelensis_HH1
  bioproject: PRJEB530
  biosample: AMD00012916
  record_type: contig
taxon:
  name: Strongyloides venezuelensis
```

(Optional) Note the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) ID for your organism (`75913`). This can be used to add extra taxonomic information to the dataset metadata.

#### Run `blobtools create`

Use the BlobTools2 `blobtools create` command to create a new [BlobDir](https://blobtoolkit.genomehubs.org/specification/) dataset in the location specified by the last argument (in this case `~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME`). The `--fasta` flag is required. The `--meta` and `--taxid` flags are optional (use of the `--taxid` flag requires `--taxdump` to provide taxonomy information).

```
~/blobtoolkit/blobtools2/blobtools create \
    --fasta ~/BTK_TUTORIAL/FILES/ASSEMBLY_NAME.fasta \
    --meta ~/BTK_TUTORIAL/FILES/ASSEMBLY_NAME.yaml \
    --taxid 75913 \
    --taxdump ~/blobtoolkit/taxdump \
    ~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME
```

This command creates a new directory with a set of files containing values for GC-content (`gc.json`), length (`length.json`), number of Ns (`ncount.json`) and sequence names (`identifiers.json`) for each sequence in the assembly. A final file (`meta.json`) contains metadata for the dataset describing the datatypes of the available fields and the ranges of values for each of these fields:

```
ls AssemblyName
gc.json identifiers.json length.json meta.json ncount.json
```

Because the `--taxid` flag was used, the metadata also includes taxonomy at eight ranks from superkingdom to species:

```
cat AssemblyName/meta.json
```

```
...
"taxon":{
  "name":"Strongyloides venezuelensis",
  "taxid":75913,
  "superkingdom":"Eukaryota",
  "kingdom":"Metazoa",
  "phylum":"Nematoda",
  "class":"Chromadorea",
  "order":"Rhabditida",
  "family":"Strongyloididae",
  "genus":"Strongyloides",
  "species":"Strongyloides venezuelensis"
 }
...
```

#### Next steps

See [_Adding data to a dataset_](https://blobtoolkit.genomehubs.org/command-line/tutorials/adding-data-to-a-dataset/) for details of how to add additional data including BLAST hits and read coverage a dataset.

See [_Updating metadata_](https://blobtoolkit.genomehubs.org/command-line/tutorials/updating-metadata/) to see how to add or update metadata after a dataset has been created.
