---
title: Adding BUSCO
date: 2020-02-12
---

Adding [BUSCO](https://busco.ezlab.org) data to a dataset allows BUSCO scores to be shown for the full dataset and recalculated for filtered datasets. Files for multiple lineages can be imported for comparison.

### Example

Example BUSCO files for _Strongyloides venezuelensis_ can be obtained from the BlobToolKit downloads site:

cd ~/BTK_TUTORIAL/FILES

curl [https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.busco.nematoda_odb9.tsv.gz](https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.busco.nematoda_odb9.tsv.gz) | gunzip -c > ASSEMBLY_NAME.busco.nematoda_odb9.tsv

curl [https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.busco.metazoa_odb9.tsv.gz](https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.busco.metazoa_odb9.tsv.gz) | gunzip -c > ASSEMBLY_NAME.busco.metazoa_odb9.tsv

These files can be imported to add BUSCO annotations to the assembly contigs:

```
~/blobtoolkit/blobtools2/blobtools add \
    --busco ~/BTK_TUTORIAL/FILES/ASSEMBLY_NAME.busco.nematoda_odb9.tsv \
    --busco ~/BTK_TUTORIAL/FILES/ASSEMBLY_NAME.busco.metazoa_odb9.tsv \
    ~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME
```

### Configuration

There are 2 configuration options when adding BUSCO results using `blobtools add`:

- `[--busco](#cov)` – BUSCO full summary file to import. _\[_**_Required_**_\]_
- `[--replace](#replace)` – Replace existing fields if present. _\[Default: false\]_

#### `--busco`

One or more BUSCO `full_summary` TSV files can be specified using the `--busco` flag. The files will be added to the dataset in the order they are specified, with the first file used to show BUSCO scores in the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/) _snail_ view.

**Usage:**

```
blobtools add \
    --busco ASSEMBLY_NAME.busco.nematoda_odb9.full_summary.tsv \
    --busco ASSEMBLY_NAME.busco.metazoa_odb9.full_summary.tsv \
    --busco ASSEMBLY_NAME.busco.eukaryota_odb9.full_summary.tsv \
    ...
```

#### `--replace`

If a `blobtools add` command would overwrite an existing field, the default behaviour is to issue a warning and not replace the existing field. To change this behaviour and allow existing fields to be overwritten, set the `--replace` flag.

**Usage:**

```
blobtools add \
    ...
    --replace \
    ...
```
