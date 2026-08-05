---
title: Getting Started with BlobTools2
date: 2019-02-11
---

BlobTools2 is a re-implementation of BlobTools designed around a more modular data structure to allow greater flexibility to incorporate additional datatypes. BlobTools2 datasets can be visualised interactively in the BlobToolKit Viewer.

To get started with BlobTools2, visit [github.com/blobtoolkit/blobtools2](https://github.com/blobtoolkit/blobtools2) for information on dependencies and to get the code.

BlobTools2 datasets are comprised of a number of files representing individual fields within a single directory (a BlobDir) containing an additional file with metadata describing the fields.

The modular nature of BlobTools2 means that a BlobDir can be created with just an assembly FASTA file. Results of analyses such as BLAST searches, read mapping and BUSCO can be added alongside the assembly FASTA, or in a subsequent command.

The guide below introduces some of the main BlobTools2 commands to create and add data to a new dataset then view that dataset interactively:

- [Set up BlobTools2](#setup_blobtools)
- [Create a BlobDir](#create_blobdir)
- [Add BLAST hits](#add_blast)
- [Add mapping coverage](#add_mapping)
- [Add BUSCO scores (optional)](#add_busco)
- [Add links (optional)](#add_links)
- [Open dataset in BlobToolKit Viewer](#open_dataset)

####

Set up BlobTools2

```
conda create -n blobtools2 -y python=3.6 docopt pyyaml ujson pysam tqdm nodejs seqtk
conda activate blobtools2
```

```
mkdir -p taxdump
cd taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
cd ..
```

```
git clone https://github.com/blobtoolkit/blobtools2
```

```
git clone https://github.com/blobtoolkit/viewer
cd viewer
npm install
cd ..
```

####

Create a BlobDir

The minimum requirement to create a new BlobDir is an assembly FASTA file:

```
./blobtools2/blobtools create \
    --fasta /path/to/assembly.fasta \
    AssemblyName
```

This command creates a new directory in the location specified by the last argument (in this case "AssemblyName") that contains a set of files containing values for GC-content (gc.json), length (length.json), number of Ns (ncount.json) and sequence names (identifiers.json) for each sequence in the assembly. A final file (meta.json) contains metadata for the dataset describing the datatypes of the available fields and the ranges of values for each of these fields:

```
ls AssemblyName
gc.json identifiers.json length.json meta.json ncount.json
```

BlobTools2 supports richer metadata that can be supplied during creation or added subsequently. Further assembly or taxonomic metadata can be added from a JSON/YAML file, for example metadata for a published assembly could include:

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

If this file is saved as "meta.yaml", it can be included in the blobtools create command as follows:

```
./blobtools2/blobtools create \
    --fasta /path/to/assembly.fasta \
    --meta meta.yaml \
    AssemblyName
```

Entries in the assembly and taxon sections of the metadata are searchable within the BlobToolKit Viewer. For taxonomic information, a shortcut to including all ranks is to specify the taxid along with the path to a local copy of the NCBI taxonomy dump:

```
./blobtools2/blobtools create \
    --fasta /path/to/assembly.fasta \
    --meta meta.yaml \
    --taxid 42157 \
    --taxdump ~/taxdump \
    AssemblyName
```

```
cat AssemblyName/meta.json
```

```
...
"taxon":{
  "genus":"Onchocerca",
  "name":"Onchocerca ochengi",
  "taxid":42157,
  "superkingdom":"Eukaryota",
  "kingdom":"Metazoa",
  "phylum":"Nematoda",
  "class":"Chromadorea",
  "order":"Rhabditida",
  "family":"Onchocercidae",
  "species":"Onchocerca ochengi"
 }
...
```

####

Add BLAST hits

The BlobTools approach uses BLAST hits to provide taxonomic annotation for each sequence in an assembly. NCBI BLAST searches must be run with an appropriate output format in order to be imported:

```
blastn -db nt \
       -query assembly.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out blast.out
```

Diamond BLAST searches should use a command similar to the following:

```
diamond blastx \
        --query assembly.fasta \
        --db /path/to/uniprot.db.with.taxids \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 16 \
        > diamond.out
```

These files can be imported to assign taxonomic labels to the assembly scaffolds using taxrules:

```
./blobtools2/blobtools add \
    --hits blast.out \
    --hits diamond.out \
    --taxrule bestsumorder \
    --taxdump ~/taxdump \
    AssemblyName
```

####

Add mapping coverage

The final data required to generate a standard blob plot is coverage information from mapping sequencing reads back to the assembly. BlobTools2 parses sorted SAM/BAM format files to calculate coverage information using the PySAM library.

To map Illumina reads using minimap2:

```
minimap2 -ax sr \
         -t 16 assembly.fasta \
         reads_1.fastq.gz reads_2.fastq.gz \
| samtools sort -@16 -O BAM -o assembly.reads.bam -
```

For PacBio or Oxford NanoPore reads, substitute `-ax sr`for `-ax map-pb`or `-ax map-ont`, respectively.

The resulting BAM file(s) can then be imported into a BlobDir:

```
./blobtools2/blobtools add \
    --cov assembly.reads.bam \
    AssemblyName
```

####

Add BUSCO scores

Incorporating BUSCO scores into a BlobDir allows them to be incorporated into dataset visualisations and to be used in assessment of the distribution of contigs/scaffolds containing BUSCO orthologues.

To run BUSCO on one or more lineages:

```
run_busco \
    -I assembly.fasta \
    -o assembly_nematoda_odb9 \
    -l /path/to/busco/lineages/nematoda_odb9 \
    -m geno \
    -c 16
```

Add the BUSCO full_table output to a BlobDir:

```
./blobtools2/blobtools add \
    --busco run_assembly_nematoda_odb9/full_table_assembly_nematoda_odb9.tsv \
    AssemblyName
```

####

Add links

Links can be added to identifiers and accessions included in the metadata, to individual records and to specific BLAST hits. These links provide access to additional resources from the detail, table and hit-distribution views in the BlobToolKit Viewer. While some of these features are more suited to public instances, adding links to BLAST hits, in particular, can be a convenient way to explore the reasons for specific taxonomic assignments in a draft assembly.

Links use a template format that allows the link to be dynamically generated based on fields in the dataset.

Add links to public assembly and taxonomy information:

```
./blobtools2/blobtools add \
    --link assembly.accession.ENA="https://www.ebi.ac.uk/ena/data/view/{accession}" \
    --link taxon.taxid.ENA="https://www.ebi.ac.uk/ena/data/view/Taxon:{taxid}" \
    --link taxon.name.Wikipedia="https://en.wikipedia.org/wiki/{name}" \
    AssemblyName
```

Add links to individual records:

```
./blobtools2/blobtools add \
    --link record.ENA="https://www.ebi.ac.uk/ena/data/view/{id}" \
    AssemblyName
```

Add links to individual BLAST hits:

```
./blobtools2/blobtools add \
    --link position.0.NCBI="https://www.ncbi.nlm.nih.gov/nuccore/{subject}" \
    --link position.1.UniProt="https://www.uniprot.org/uniprot/{subject}" \
    AssemblyName
```

####

Open dataset in BlobToolKit Viewer

The BlobToolKit Viewer API loads all BlobDir contained within a given directory and makes the data available to the Viewer. To start the API and viewer to host all datasets in the current directory:

```
./blobtools2/blobtools host `pwd`
```

```
Starting BlobToolKit API on port 8081 (pid: 1234)
Starting BlobToolKit viewer on port 8080 (pid: 1235)
Visit http://localhost:8080 to use the interactive BlobToolKit Viewer.
```

Visit the
