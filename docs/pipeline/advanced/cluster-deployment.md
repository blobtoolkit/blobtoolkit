---
layout: default
title: Cluster Deployment
---

# Running the Pipeline on a Cluster

## Preferred: Nextflow Pipeline

The recommended way to run BlobToolKit on a cluster is via the
[sanger-tol/blobtoolkit](https://github.com/sanger-tol/blobtoolkit) Nextflow
pipeline. Nextflow handles job submission, resource management and
containerisation natively.

See the pipeline documentation at
https://pipelines.tol.sanger.ac.uk/blobtoolkit for full parameter details and
profile options.

### Example: LSF cluster with Singularity

```bash
nextflow run sanger-tol/blobtoolkit \
    -r 0.10.0 \
    -profile sanger,singularity \
    --input samplesheet.csv \
    --fasta assembly.fasta \
    --accession GCA_000000000.1 \
    --taxon 9606 \
    --outdir results/ \
    --taxdump /data/taxonomy/new_taxdump \
    --blastp /data/uniprot/reference_proteomes.dmnd \
    --blastn /data/nt/nt.nal \
    --blastx /data/uniprot/reference_proteomes.dmnd \
    --align
```

Or using the `btk pipeline run` wrapper once your config file is set up
(see [Configuration](../configuration)):

```bash
btk pipeline run --nextflow --config pipeline.yaml
```

### Selecting a Profile

Nextflow profiles bundle executor, container engine and resource settings.
Common options:

| Profile                | Description                                    |
| ---------------------- | ---------------------------------------------- |
| `sanger,singularity`   | Sanger Institute LSF cluster with Singularity. |
| `standard,singularity` | Local executor with Singularity.               |
| `standard,docker`      | Local executor with Docker.                    |

See the pipeline docs for a full list of available profiles and how to write a
custom institutional profile.

### Resuming a Run

Nextflow automatically caches completed tasks. To resume after a failure:

```bash
nextflow run sanger-tol/blobtoolkit -r 0.10.0 -resume ...
```

---

## Legacy: Snakemake Pipeline

> **Note**: The Snakemake-based pipeline is no longer actively maintained.
> It is documented here for users running existing workflows.

The legacy pipeline was a Snakemake workflow submitted via a job script. Key
variables to set in your submission script:

```bash
export PIPELINE=/path/to/blobtoolkit/insdc-pipeline
export WORKDIR=/path/to/workdir          # must contain $ASSEMBLY.yaml
export CONDA_DIR=/path/to/.conda
export THREADS=32
```

Run Snakemake inside a job script:

```bash
eval "$(conda shell.bash hook)"
conda activate btk_env

snakemake -p \
    --use-conda \
    --conda-prefix "$CONDA_DIR" \
    --directory "$WORKDIR/" \
    --configfile "$WORKDIR/$ASSEMBLY.yaml" \
    --latency-wait 60 \
    --rerun-incomplete \
    -j "$THREADS" \
    -s "$PIPELINE/Snakefile" \
    --resources btk=1
```

### Submitting Per-Rule Jobs (drmaa)

On clusters with drmaa support, Snakemake can submit each rule as a separate
job using a `cluster.yaml` resource file:

```bash
snakemake -p \
    --use-conda \
    --conda-prefix "$CONDA_DIR" \
    --directory "$WORKDIR/" \
    --configfile "$WORKDIR/$ASSEMBLY.yaml" \
    --cluster-config "$CLUSTER_CONFIG" \
    --drmaa " -o {log}.o -e {log}.e \
              -R 'select[mem>{cluster.mem}] rusage[mem={cluster.mem}]' \
              -M {cluster.mem} -n {cluster.threads} -q {cluster.queue}" \
    --latency-wait 60 --rerun-incomplete \
    -j "$THREADS" \
    -s "$PIPELINE/Snakefile" \
    --resources btk=1
```

An example `cluster.yaml` is available in the
[insdc-pipeline repository](https://github.com/blobtoolkit/insdc-pipeline).
