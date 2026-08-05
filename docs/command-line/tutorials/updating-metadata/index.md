---
title: Updating metadata
date: 2020-02-11
---

Dataset metadata, including links to external resources, can be added or updated using the `blobtools add` command. Given the focus on updating metadata, all examples use `blobtools replace` which is equivalent to `blobtools add --replace`.

#### Updating assembly/taxon information

Any values in the BlobDir dataset metadata can be updated by specifying `--key path=value`. In this context, `path` is a `.`\-separated hierarchy of keys, e.g. `assembly.accession`.

~/blobtoolkit/blobtools2/blobtools replace \\
--key assembly.reference=doi:10.1038/ng.3495 \\
--key taxon.common_name=threadworm \\
~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME

#### Updating default plot axes

The default plot axes when a dataset is opened in the BlobToolKit Viewer, are defined in the `plot` section of the metadata. This will usually be:

- x – GC proportion
- y – First coverage field added (see [Adding coverage](https://blobtoolkit.genomehubs.org/command-line/tutorials/adding-data-to-a-dataset/adding-coverage/))
- z – Length
- cat (categories) – Taxonomy from the first taxrule added at the phylum level (see [Adding hits](https://blobtoolkit.genomehubs.org/command-line/tutorials/adding-data-to-a-dataset/adding-hits/))

To update this, use the `--key` option as above:

~/blobtoolkit/blobtools2/blobtools replace \\
--key plot.x=length \\
--key plot.y=DRR008460_cov \\
--key plot.cat=bestsumorder_family \\
~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME

#### Adding links to external resources

Datasets can contain links to external resources, these are shown in the _detail_ and _table_ views of the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/). Links are defined using URL templates that allow dataset specific terms to be included in the URL. The option to add new links (`--link path=URL_template`) uses `path` in the same way as for `--key` above for values in the metadata (e.g. `assembly.accession`), with additional values to allow links from individual contains or BLAST hits. A set of default links is available when starting the Viewer. When specified using JSON (in a file named `default.json` in the same directory as the datasets) these links can have custom titles. When added to the dataset, however, the link title is set using the last key in the `path`:

\--link 'taxon.taxid=https://www.ebi.ac.uk/ena/browser/view/Taxon:{taxid}'
--link 'taxon.species=https://wikipedia.org/{species}'
--link 'assembly.bioproject=https://www.ebi.ac.uk/ena/browser/view/{bioproject}'

For comparison, more complex linking relationships can be specified in `default.json` or by directly editing the dataset `meta.json` file.

{
...
"links": {
"taxon": {
"taxid" : {
"ENA": "https://www.ebi.ac.uk/ena/browser/view/Taxon:{taxid}"
},
"species" : {
"Wikipedia": "https://wikipedia.org/{species}"
}
},
"blobtoolkit": {
"commit" : {
"Github": "{pipeline}/tree/{commit}"
}
},
"assembly": {
"biosample": {
"ENA": "https://www.ebi.ac.uk/ena/browser/view/{biosample}"
},
"bioproject": {
"ENA": "https://www.ebi.ac.uk/ena/browser/view/{bioproject}"
}
},
"position": \[
{
"patterns": \[
{
"title": "NCBI RefSeq",
"template": "https://www.ncbi.nlm.nih.gov/nuccore/{subject}",
"regex": "^\[NXW\]\[A-Z\]\_.+"
},
{
"title": "ENA",
"template": "https://www.ebi.ac.uk/ena/browser/view/{subject}",
"regex": "^.+$"
}
\]
},
{
"UniProt": "https://www.uniprot.org/uniprot/{subject}"
}
\]
}
...
}

For links based on information in the metadata, the link will only be added if the link URL can be resolved. Adding the `--skip-link-test` flag will allow links to be added even if the URL cannot be resolved.
