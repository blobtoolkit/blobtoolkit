---
title: Adding coverage
date: 2020-02-12
---

After loading `--hits`, the final data required to generate a standard blob plot is coverage information from mapping sequencing reads back to the assembly. BlobTools2 parses sorted BAM/SAM/CRAM format files to calculate coverage information using the PySAM library.

### Example

Due to the size of read alignment files, BAM files for analysed datasets on the [public Viewer instance](https://blobtoolkit.genomehubs.org/view) are not retained. To add coverage information for _Strongyloides venezuelensis_ without recomputing the alignment, it is possible to import coverage from a [BlobDir](https://blobtoolkit.genomehubs.org/specification/) dataset. While this is not the typical usage, the command remains the same as the filetype (in this case JSON) is inferred from the file extension.

Fetch and extract coverage data for _Strongyloides venezuelensis_ from the BlobToolKit downloads site:

curl -L [https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.blobdir.tar.gz](https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.blobdir.tar.gz) | tar -xJf - S_venezuelensis_HH1/DRR008460_cov.json

This file can be imported to add coverage information to the assembly contigs:

~/blobtoolkit/blobtools2/blobtools add \\  
 --cov ~/BTK_TUTORIAL/FILES/S_venezuelensis_HH1/DRR008460_cov.json \\  
 ~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME

### Configuration

There are 3 configuration options when adding coverage using `blobtools add`:

- `[--cov](#cov)` – BAM/SAM/CRAM read alignment file to import. _\[_**_Required_**_\]_
- `[--threads](#threads)` – Number of threads to use when parsing read alignment files. _\[Default: 1\]_
- `[--replace](#replace)` – Replace existing fields if present. _\[Default: false\]_

#### `--cov`

In the more typical usage, one or more BAM, SAM or CRAM format read mapping files can be specified using the `--cov` flag. The [BlobToolKit Pipeline](https://blobtoolkit.genomehubs.org/pipeline/) uses [minimap2](https://github.com/lh3/minimap2) to generate these files, with commands similar to:

```
minimap2 -ax sr \
         -t 16 ASSEMBLY_NAME.fasta \
         DRR008460_1.fastq.gz DRR008460_2.fastq.gz \
| samtools sort -@16 -O BAM -o ASSEMBLY_NAME.DRR008460.bam -
```

For PacBio or Oxford NanoPore reads, substitute `-ax sr`for `-ax map-pb`or `-ax map-ont`, respectively.

When importing, the file basename is used to determine the resulting field names so `ASSEMBLY_NAME.DRR008460.bam` would generate a base coverage field named `ASSEMBLY_NAME.DRR008460_cov` and a read coverage field named `ASSEMBLY_NAME.DRR008460_read_cov`. To set a different name, use the syntax `--cov FILE_NAME.bam=FIELD_NAME` to add `FIELD_NAME_cov` and `FIELD_NAME_read_cov` fields.

**Usage:**

```
blobtools add \
    --cov ASSEMBLY_NAME.DRR008460.bam \
    --cov ASSEMBLY_NAME.DRR008461.bam=DRR008461 \
    ...
```

#### `--threads`

BlobTools2 can use multiple threads when parsing read coverage files. To use more than one thread, set `--threads` to an appropriate value.

**Usage:**

```
blobtools add \
    ...
    --threads 8 \
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
