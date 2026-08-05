---
title: Adding data to a dataset
date: 2020-02-10
---

Once you have a [BlobDir](https://blobtoolkit.genomehubs.org/specification/) dataset (see [Creating a dataset](https://blobtoolkit.genomehubs.org/command-line/tutorials/creating-a-dataset/)), further data can be added by parsing analysis output files into one or more fields using the `blobtools add` command. Dedicated parsers are available for a range of analysis types, assigning values to each contig in an assembly: BLAST or Diamond hits provide taxonomic assignments; read mapping files provide base and read coverage; BUSCO results show completeness metrics for full and filtered assemblies; and fields from generic text files can be imported for maximum flexibility.

All `blobtools add` options can be passed directly to `blobtools create` to allow dataset creation and analysis import in a single step.

#### Adding hits

The BlobTools approach uses BLAST hits to provide taxonomic annotation for each sequence in an assembly. When run using the BlobToolKit Pipeline, a wrapper script is used to break long sequences into chunks to obtain a distribution of BLAST hits and closely related taxa can be automatically filtered out. This is filtering step is particularly important for publicly available assemblies as the BLAST databases may already contain the query sequence. [Read more...](https://blobtoolkit.genomehubs.org/command-line/tutorials/adding-data-to-a-dataset/adding-hits/)

#### Adding coverage

After loading `--hits`, the final data required to generate a standard blob plot is coverage information from mapping sequencing reads back to the assembly. BlobTools2 parses sorted BAM/SAM/CRAM format files to calculate coverage information using the PySAM library. [Read more...](https://blobtoolkit.genomehubs.org/adding-coverage/)

#### Adding BUSCO

Adding [BUSCO](https://busco.ezlab.org) data to a dataset allows BUSCO scores to be shown for the full dataset and recalculated for filtered datasets. Files for multiple lineages can be imported for comparison. [Read more...](https://blobtoolkit.genomehubs.org/command-line/tutorials/adding-data-to-a-dataset/adding-busco/)

#### Adding text files

Files can be imported from generic text files, with the ability to specify column separators and map columns to field names. This provides flexibility to view and filter a wide range of additional analysis outputs beyond those with dedicated parsers. [Read more...](https://blobtoolkit.genomehubs.org/command-line/tutorials/adding-data-to-a-dataset/adding-text-files/)
