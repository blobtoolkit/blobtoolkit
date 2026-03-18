---
layout: default
title: Pipeline Configuration
---

# Pipeline Configuration

The preferred BlobToolKit pipeline configuration is documented externally at:

- https://pipelines.tol.sanger.ac.uk/blobtoolkit

This page is retained as a legacy/local configuration reference for existing scripts and older workflows.

Some fields below were introduced for Snakemake-based workflows and may not be required for modern deployments.

## Configuration File Structure

Typical configuration files are divided into these sections:

- `assembly`
- `busco`
- `reads`
- `settings`
- `similarity`
- `taxon`
- `keep_intermediates`

Example:

```yaml
assembly:
	accession: GCA_000298335.1
	alias: DroAlb_1.0
	level: scaffold
	scaffold-count: 26354
	span: 253560284
	prefix: ACVV01

busco:
	lineage_dir: /path/to/busco/lineages
	lineages:
		- diptera_odb10
		- arthropoda_odb10
		- eukaryota_odb10

reads:
	paired:
		-
			- SRR026696
			- ILLUMINA
			- 482114248
			- ftp.sra.ebi.ac.uk/vol1/fastq/SRR026/SRR026696/SRR026696_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR026/SRR026696/SRR026696_2.fastq.gz
	single: []
	coverage:
		max: 100
		min: 0.5

settings:
	taxonomy: /path/to/taxdump
	tmp: /tmp
	blast_chunk: 100000
	blast_max_chunks: 10
	blast_overlap: 500

similarity:
	defaults:
		evalue: 1e-25
		max_target_seqs: 10
		root: 1
		mask_ids: []
	databases:
		-
			local: /path/to/ncbi_nt
			name: nt
			source: ncbi
			tool: blast
			type: nucl
		-
			local: /path/to/uniprot
			name: reference_proteomes
			source: uniprot
			tool: diamond
			type: prot
	taxrule: bestsumorder

taxon:
	taxid: 7291
	name: Drosophila albomicans

keep_intermediates: true
```

## Section Notes

### assembly

`assembly` defines core assembly metadata, accession context, and file naming prefix.

- `prefix` is typically used as the assembly identifier across outputs.
- For local/non-public assemblies, `accession` is often set to `draft`.
- `span` can be used to estimate and control coverage-aware read subsampling.

### busco

`busco` controls BUSCO lineages and lineage database location.

- Set `lineages: []` to skip BUSCO runs.
- Use one or more lineage datasets depending on expected taxonomic placement.

### reads

`reads` declares paired and/or single-end libraries.

- Paired and single sections can both be used.
- Platforms such as `ILLUMINA`, `OXFORD_NANOPORE`, and `PACBIO_SMRT` affect mapping defaults.
- Optional `coverage.max` and `coverage.min` values can be used to gate or subsample read processing.

### settings

`settings` includes runtime and search tuning options.

- `taxonomy` should point to a compatible NCBI taxdump.
- `tmp` should reference a location with sufficient disk space.
- Chunk-related fields (`blast_chunk`, `blast_max_chunks`, `blast_overlap`) tune long-sequence BLAST partitioning.

### similarity

`similarity` configures sequence search behavior and databases.

- `defaults` apply globally.
- `mask_ids` can exclude clades (for example, to avoid self-hits in public assemblies).
- `databases` lists search backends and local paths.
- `taxrule` is usually `bestsumorder`.

### taxon

`taxon` should include both the NCBI taxonomy ID and display name for metadata and downstream categorisation.

### keep_intermediates

Set `keep_intermediates: true` to retain temporary/intermediate files for debugging and reuse.

## Legacy Snakemake Context

Older guides reference direct Snakemake execution with environment-specific variables. That material is still useful in advanced setups, but preferred production usage should follow https://pipelines.tol.sanger.ac.uk/blobtoolkit.

## Running With a Configuration

Preferred run/config route:

- https://pipelines.tol.sanger.ac.uk/blobtoolkit

Legacy/local CLI route:

```bash
btk pipeline run --help
```

## Database Setup

Database preparation is documented in [Installation Databases](../installation/databases).

## Samplesheet Format (If Required)

Some workflow variants may still expect tabular sample declarations.

Define your samples in a samplesheet file:

```
sample,reads1,reads2,type
sample1,reads1_R1.fq,reads1_R2.fq,paired
sample2,reads2_R1.fq,reads2_R2.fq,paired
sample3,reads3.fq,,single
```

## Configuration Options

Refer to https://pipelines.tol.sanger.ac.uk/blobtoolkit for current options and defaults.

For legacy/local workflows, use btk pipeline help and your installed release notes.
