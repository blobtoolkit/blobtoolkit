---
layout: default
title: Database Setup
---

# Database Setup

BlobToolKit workflows that use taxonomy and sequence similarity require local reference databases.

This page provides practical shell commands for setting up commonly used databases.

## 1. NCBI Taxdump

```bash
mkdir -p taxdump
cd taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
cd -
```

## 2. NCBI nt

```bash
mkdir -p nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz" -P nt/
for file in nt/*.tar.gz; do
  tar xf "$file" -C nt && rm "$file"
done
```

## 3. UniProt Reference Proteomes

The UniProt setup below builds a Diamond database and taxon map.

```bash
mkdir -p uniprot
wget -q -O uniprot/reference_proteomes.tar.gz \
  ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
    -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
    awk '/tar.gz/ {print $9}')

cd uniprot
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

printf "accession\taccession.version\ttaxid\tgi\n" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb \
  -p 16 \
  --in reference_proteomes.fasta.gz \
  --taxonmap reference_proteomes.taxid_map \
  --taxonnodes ../taxdump/nodes.dmp \
  --taxonnames ../taxdump/names.dmp \
  -d reference_proteomes.dmnd
cd -
```

## 4. BUSCO Lineages

Download the lineages required for your analyses.

Example:

```bash
mkdir -p busco
wget -q -O eukaryota_odb10.tar.gz "https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz"
tar xf eukaryota_odb10.tar.gz -C busco
```

## Notes

- Keep taxdump and sequence databases version-aligned where possible.
- Store databases in stable paths and reference those paths in pipeline configuration.
- Database creation can require substantial CPU, RAM, and disk space.
