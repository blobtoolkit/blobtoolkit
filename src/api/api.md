# API Documentation

Data and pre-generated static images can be retrieved directly via the BlobToolKit viewer API.

## Search

### /api/v1/search/all

Returns an array of top-level dataset metadata containing dataset `name`, `id` and any additional searchable fields for all hosted datasets.

```
curl "https://example.com/api/v1/search/all"
```

```json
[
  {
    "id":"AAGD02",
    "name":"AAGD02",
    "accession":"GCA_000149515.1",
    "species":"",
    ...},
  {
    "id":"ADBU02",
    "name":"ADBU02",
    ...
  },
  ...
]
```

### /api/v1/search/:term

Returns an array of top-level dataset metadata for all datasets matching `:term`.

```
curl "https://example.com/api/v1/search/Caenorhabditis"
```

```json
[
  {
    "genus":"Caenorhabditis",
    "id":"AAGD02",
    ...},
  {
    "genus":"Caenorhabditis",
    "id":"LFJK02",
    ...
  },
  ...
]
```

## Dataset

### /api/v1/dataset/id/:dataset_id

Returns full metadata for the dataset with the provided `:dataset_id`.

```
curl "https://example.com/api/v1/dataset/id/LFJK02"
```

```json
{
  "id":"LFJK02",
  "assembly":{
    "accession":"GCA_001643735.2",
    "bioproject":"PRJNA248909",
    "level":"scaffold",
    "scaffold-count":1591,
    "span":118549266,
    ...
  },
  "fields":[
    {
      "id":"identifiers",
      "type":"identifier"
    },
    {
      "id":"gc",
      "datatype":"float",
      "range":[0.2498,0.6325],
      ...
    },
    {
      "id":"length",
      ...
    },
    ...
  ],
  "reads":{
    "SRR5837619":{
      "Mapped reads":"52760532",
      "Total reads":"59124810",
      "coverage":502.1056381572,
      ...
    },
    "SRR5837620":{
      ...
    },
    ...
  },
  "record_type":"scaffold",
  "records":1591,
  "settings":{
    "commit":"843f7af",
    "pipeline":"https://github.com/blobtoolkit/insdc-pipeline.git",
    ...
  },
  "static_plots":"true",
  "taxon":{
    "genus":"Caenorhabditis",
    "taxid":31234,
    "phylum":"Nematoda",
    "species":"Caenorhabditis remanei",
    ...
  }
}
```

### /api/v1/dataset/id/:dataset_id/:key

Returns the portion of the dataset metadata matching `:key`.

```
curl "https://example.com/api/v1/dataset/id/LFJK02/assembly"
```

```json
{
  "accession":"GCA_001643735.2",
  "alias":"ASM164373v2",
  "bioproject":"PRJNA248909",
  "biosample":"SAMN02803904",
  "level":"scaffold",
  "prefix":"LFJK02",
  "scaffold-count":1591,
  "span":118549266,
  "file":"LFJK02.fasta"
}
```

```
curl "https://example.com/api/v1/dataset/id/LFJK02/static_plots"
```

```json
"true"
```

## Summary

### /api/v1/summary/:dataset_id

Returns a pre-computed summary for `:dataset_id`, if available.

```
curl "https://example.com/api/v1/summary/AACW02"
```

```json
{
  "id": "current",
  "datasetId": "AACW02",
  "taxon": "Rhizopus delemar RA 99-880",
  "taxid": 246409,
  "search": "all",
  "params": {},
  "summaryStats": {
    "hits": {
      "total": {
        "span": 46148878,
        "count": 83,
        "n50": 3104119,
        "l50": 6,
        "n90": 783160,
        "l90": 17
      },
      "Mucoromycota": {
        "count": 68,
        "span": 42389618,
        "n50": 2712543
      },
      "no-hit": {
        "count": 4,
        "span": 28432,
        "n50": 7630
      },
      ...
    },
    "taxonomy": {
      "taxid": 246409,
      "lineage": "Eukaryota; Fungi; Mucoromycota; Mucoromycetes; Mucorales; Rhizopodaceae; Rhizopus; Rhizopus delemar",
      "target": "Mucoromycota",
      "targetRank": "phylum"
    },
    "baseComposition": {
      "at": 0.644,
      "gc": 0.356,
      "n": 0.0178
    },
    "busco": {
      "fungi_odb9": {
        "c": 276,
        "d": 112,
        "m": 9,
        "f": 5,
        "t": 290,
        "s": 164,
        "string": "C:95.2%[S:56.6%,D:38.6%],F:1.7%,M:3.1%,n:290"
      },
      ...
    },
    "stats": {
      "noHit": 0.000616,
      "target": 0.919,
      "spanOverN50": 14.9
    }
  }
}
```

## Image

The following routes are available if static images have been pre-generated for a dataset (i.e. if `/dataset/id/:dataset_id/static_plots` returns `"true"`)

### /api/v1/image/:dataset_id

Returns the default image for a dataset, which is a `PNG` format square-binned blobplot at maximum available resolution.

```
curl -o LFJK02.png "https://example.com/api/v1/image/LFJK02"
```

### /api/v1/image/:dataset_id/:view

Returns an image for a specified dataset view.

```
curl -o LFJK02.blob.png "https://example.com/api/v1/image/LFJK02/blob"
```
```
curl -o LFJK02.cumulative.png "https://example.com/api/v1/image/LFJK02/cumulative"
```
```
curl -o LFJK02.snail.png "https://example.com/api/v1/image/LFJK02/snail"
```

### /api/v1/image/:dataset_id/:view/:type

Returns a specific `:type` for dataset views with more then one static image type.

```
curl -o LFJK02.blob.circle.png "https://example.com/api/v1/image/LFJK02/blob/circle"
```
```
curl -o LFJK02.blob.hex.png "https://example.com/api/v1/image/LFJK02/blob/hex"
```
```
curl -o LFJK02.blob.square.png "https://example.com/api/v1/image/LFJK02/blob/square"
```

### query string parameters

#### /api/v1/image/:dataset_id/:type?format=:format

Determines the image `:format` to retrieve, typically `PNG` and `SVG` formats will be available. If not specified, the default format is `PNG`.

```
curl -o LFJK02.blob.png "https://example.com/api/v1/image/LFJK02/blob?format=png"
```
```
curl -o LFJK02.blob.svg "https://example.com/api/v1/image/LFJK02/blob?format=svg"
```

#### /api/v1/image/:dataset_id/:type?width=:width

For `PNG` format images, determines the `:width` (in pixels) of the image to be returned. If the width specified is greater than the width of the available image, the image will be returned at its native resolution and will not be scaled up.

```
curl -o LFJK02.cumulative.1000.svg "https://example.com/api/v1/image/LFJK02/cumulative?width=1000"
```
