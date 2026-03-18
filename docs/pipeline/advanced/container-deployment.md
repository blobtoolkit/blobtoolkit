---
layout: default
title: Container Deployment
---

# Running the Pipeline in a Container

## Preferred: Nextflow with Singularity or Docker

The [sanger-tol/blobtoolkit](https://github.com/sanger-tol/blobtoolkit)
Nextflow pipeline manages container execution automatically via profiles.
Pass `-profile singularity` or `-profile docker` to select the container
engine — no manual `docker run` or `singularity exec` calls are needed.

```bash
nextflow run sanger-tol/blobtoolkit \
    -r 0.10.0 \
    -profile standard,singularity \
    --input samplesheet.csv \
    --fasta assembly.fasta \
    --accession GCA_000000000.1 \
    --taxon 9606 \
    --outdir results/
```

See [Cluster Deployment](cluster-deployment) for a full example with all
database flags and the `btk pipeline run` wrapper.

---

## Legacy: Snakemake Pipeline inside Docker

> **Note**: The Snakemake/Docker approach is no longer actively maintained.
> It is documented here for reference.

When running the full legacy pipeline inside Docker, all file paths in the
configuration file must be **relative to the container filesystem** rather than
your local filesystem, because Snakemake itself runs inside the container.

This differs from the cluster approach where Snakemake runs outside the
container and maps paths automatically.

### Example configuration (container-relative paths)

```yaml
assembly:
  accession: draft
  level: scaffold
  prefix: MyAssembly
busco:
  lineages:
    - diptera_odb10
  lineage_dir: /blobtoolkit/databases/busco
reads:
  paired:
    - - illumina_reads
      - ILLUMINA
settings:
  taxonomy: /blobtoolkit/databases/ncbi_taxdump
  tmp: /tmp
similarity:
  databases:
    - local: /blobtoolkit/databases/ncbi_db
      name: nt
      source: ncbi
      tool: blast
      type: nucl
    - local: /blobtoolkit/databases/uniprot_db
      name: reference_proteomes
      source: uniprot
      tool: diamond
      type: prot
taxon:
  taxid: 7291
  name: Drosophila albomicans
```

### Running the container

```bash
WORKDIR=/path/to/working/directory
BUSCO_DIR=/path/to/busco_lineages
UNIPROT_DIR=/path/to/uniprot
NT_DIR=/path/to/ncbi_nt
TAXDUMP_DIR=/path/to/ncbi_taxdump
CONDA_DIR=/path/to/.conda
ASSEMBLY=MyAssembly
THREADS=32

docker run -it --rm \
    -u $UID:$GROUPS \
    -v $WORKDIR:/blobtoolkit/data \
    -v $BUSCO_DIR:/blobtoolkit/databases/busco \
    -v $UNIPROT_DIR:/blobtoolkit/databases/uniprot_db \
    -v $NT_DIR:/blobtoolkit/databases/ncbi_db \
    -v $TAXDUMP_DIR:/blobtoolkit/databases/ncbi_taxdump \
    -v $CONDA_DIR:/blobtoolkit/.conda \
    genomehubs/blobtoolkit:latest \
    snakemake -p \
        --use-conda \
        --conda-prefix /blobtoolkit/.conda \
        --directory /blobtoolkit/data \
        --configfile /blobtoolkit/data/$ASSEMBLY.yaml \
        -j $THREADS \
        -s /blobtoolkit/insdc-pipeline/Snakefile \
        --resources btk=1
```

### Key points

- Volume mounts (`-v`) map local directories to the paths used in your
  configuration YAML.
- Database directories must be mounted and referenced at their
  **container-side** paths in the config.
- `-u $UID:$GROUPS` ensures output files are owned by your local user.
- The `genomehubs/blobtoolkit` image is available on
  [Docker Hub](https://hub.docker.com/r/genomehubs/blobtoolkit).
