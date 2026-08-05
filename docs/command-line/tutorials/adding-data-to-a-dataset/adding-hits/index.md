---
title: Adding hits
date: 2020-02-12
---

The BlobTools approach uses BLAST hits to provide taxonomic annotation for each sequence in an assembly. When run using the BlobToolKit Pipeline, a wrapper script is used to break long sequences into chunks to obtain a distribution of BLAST hits and closely related taxa can be automatically filtered out. This is particularly important for publicly available assemblies as the BLAST databases may already contain the query sequence.

### Example

Example BLAST and Diamond BLAST hit files for _Strongyloides venezuelensis_ can be obtained from the BlobToolKit downloads site:

cd ~/BTK_TUTORIAL/FILES

curl [https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.blastn.nt_v5.root.1.minus.6247.out.gz](https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.blastn.nt_v5.root.1.minus.6247.out.gz) | gunzip -c > ASSEMBLY_NAME.ncbi.blastn.out

curl [https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.diamond.reference_proteomes.root.1.minus.6247.out.gz](https://blobtoolkit.genomehubs.org/download/ALIAS/S_v/S_venezuelensis_HH1/S_venezuelensis_HH1.diamond.reference_proteomes.root.1.minus.6247.out.gz) | gunzip -c > ASSEMBLY_NAME.diamond.blastx.out

These files can be imported to assign taxonomic labels to the assembly contigs:

```
~/blobtoolkit/blobtools2/blobtools add \
    --hits ~/BTK_TUTORIAL/FILES/ASSEMBLY_NAME.ncbi.blastn.out \
    --hits ~/BTK_TUTORIAL/FILES/ASSEMBLY_NAME.diamond.blastx.out \
    --taxrule bestsumorder \
    --taxdump ~/blobtoolkit/taxdump \
    ~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME
```

This example uses the `--taxrule bestsumorder`, which assigns taxonomy based on the sum of bitscores for the hits to a given taxon based on the first file, only using subsequent results files for contigs with no hits a previous file.

### Configuration

There are a number of configuration options when adding hits using `blobtools add`:

- `[--hits](#hits)` – BLAST or Diamond blast file to import. _\[_**_Required_**_\]_
- `[--hits-cols](#hits-cols)` – Specify the BLAST file column order. _\[Default: 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=sstart,11=send,14=evalue\]_
- `[--taxrule](#taxrule)` – Rule to use when assigning BLAST hits to taxa. _\[Default: bestsumorder\]_
- `[--taxdump](#taxdump)` – Directory containing NCBI taxdump files. _\[_**_Required_**_\]_
- `[--evalue](#evalue)` – Set evalue cutoff when parsing hits file. _\[Default: 1\]_
- `[--bitscore](#bitscore)` – Set bitscore cutoff when parsing hits file. _\[Default: 1\]_
- `[--hit-count](#hit-count)` – Number of hits to parse when inferring taxonomy. _\[Default: 1\]_
- `[--replace](#replace)` – Replace existing fields if present. _\[Default: false\]_

#### `--hits`

One or more BLAST or Diamond BLAST tabular output format files can be specified using the `--hits` flag. The order of columns in these files can be defined using `[--hits-cols](#hits-cols)` (described below), but the default import format can be obtained with commands similar to:

**NCBI blastn against the nt database** (see [Install BlobToolKit > Databases](https://blobtoolkit.genomehubs.org/install/#databases)):

```
blastn -db nt \
       -query ASSEMBLY_NAME.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out ASSEMBLY_NAME.ncbi.blastn.out"
```

**Diamond blast against UniProt** (see [Install BlobToolKit > Databases](https://blobtoolkit.genomehubs.org/install/#databases)):

```
diamond blastx \
        --query ASSEMBLY_NAME.fasta \
        --db /path/to/uniprot.db.with.taxids \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 16 \
        > ASSEMBLY_NAME.diamond.blastx.out
```

**Usage:**

```
blobtools add \
    --hits ASSEMBLY_NAME.ncbi.blastn.out \
    --hits ASSEMBLY_NAME.diamond.blastx.out \
    ...
```

#### `--hits-cols`

Tablular BLAST output files with a different column order can be imported by specifying a comma separated list of `<column_number=field_name>`. The default value of `1=qseqid,2=staxids,3=bitscore,5=sseqid,10=qstart,11=qend,14=evalue` matches the output from the example commands above. Columns not listed in the default setting are not used.

**Usage:**

```
blobtools add \
    ...
    --hits-cols 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=qstart,11=qend,14=evalue
    ...
```

#### `--taxrule`

BlobTools2 assigns a putative taxonomic origin to each scaffold at 8 ranks from superkingdom to species based on aggregating hits in the `[--hits](#hits)` files. The `--taxrule` flag determines the rule to use when assigning BLAST hits to taxa. Two options are available: `bestsum` and `bestsumorder`. Of the taxa represented in the BLAST hits, `bestsum` assigns the taxon for which the sum of bitscores is greatest across all files while the default `bestsumorder` taxrule uses the maximum bitscore in the first file and only uses hits from subsequent files for taxa without a hit in previous files.

**Usage:**

```
blobtools add \
    ...
    --taxrule bestsum \
    ...
```

#### `--taxdump`

Taxonomy information must be loaded from a local copy of the NCBI taxdump (see [Install BlobToolKit > Databases](https://blobtoolkit.genomehubs.org/install/#databases)).

**Usage:**

```
blobtools add \
    ...
    --taxdump ~/blobtoolkit/taxdump \
    ...
```

#### `--evalue`

BLAST results can be filtered based on an evalue cutoff that will be applied before the `[--taxrule](#taxrule)`. Any hits with an evalue weaker than the value specified will be excluded.

**Usage:**

```
blobtools add \
    ...
    --evalue 1e-50 \
    ...
```

#### `--bitscore`

BLAST results can be filtered based on a bitscore cutoff that will be applied before the `[--taxrule](#taxrule)`. Any hits with an bitscore lower the value specified will be excluded.

**Usage:**

```
blobtools add \
    ...
    --bitscore 100 \
    ...
```

#### `--hit-count`

By default the 10 highest scoring hits to a given taxon will be used when applying the `[--taxrule](#taxrule)`. To use more or fewer hits, set the `--hit-count` accordingly.

**Usage:**

```
blobtools add \
    ...
    --hit-count 5 \
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
