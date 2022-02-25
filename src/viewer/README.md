# BlobToolKit Viewer (v3.0.0)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1134794.svg)](https://doi.org/10.5281/zenodo.1134794)

## Upgrading from v2.x.x

From v3.0.0, the viewer and api are available as packaged executables. These will be compatible with other v3 BlobToolKit repositories as they are released. Meanwhile, please use the latest `release/v2.6.x` branch as part of the pipeline or with `blobtools view`. To run the v3 viewer, fetch the executables from the release page and add them to your PATH.

To view a demo BlobDir, run the following commands in separate terminal windows:

```
blobtoolkit-api
```

```
blobtoolkit-viewer
```

Then open [localhost:8080/view/all](http://localhost:8080/view/all) in a web browser.

Port numbers and data directories can be set using environment variables. To view BlobDirs in a local directory and change the port numbers, run:

```
BTK_API_PORT=8880 BTK_PORT=8881 BTK_FILE_PATH=/Users/rchallis/projects/blobtoolkit/datasets blobtoolkit-api

BTK_API_PORT=8880 BTK_PORT=8881 blobtoolkit-viewer
```

## About BlobToolKit Viewer

BlobToolKit Viewer is a genome-scale dataset visualistion tool developed as part of the [blobtoolkit](https://blobtoolkit.genomehubs.org) project to allow browser-based identification and filtering of target and non-target data in genome assemblies.

We are running the BlobToolKit [insdc-pipeline](https://github.com/blobtoolkit/insdc-pipeline) on all public (INSDC registered) eukaryote genome assemblies and making the results available in an instance of this viewer at [blobtoolkit.genomehubs.org](https://blobtoolkit.genomehubs.org/view/).
