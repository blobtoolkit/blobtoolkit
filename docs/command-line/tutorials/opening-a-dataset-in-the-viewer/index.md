---
title: Opening a dataset in the Viewer
date: 2020-02-10
---

Datasets created using [BlobTools2](https://blobtoolkit.genomehubs.org/blobtools2/) can be viewed interactively by starting a local instance of the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/btk-viewer/). The Viewer can be started using the `blobtools host` and `blobtools view` commands. Of these, `blobtools view` offers more flexibility and can also be used to generate plots non-interactively (see _[Generating plots on the command line](https://blobtoolkit.genomehubs.org/command-line/tutorials/generating-plots-on-the-command-line/)_).

### Example

Host a dataset on your local computer or a remote server and ready to view in a browser on your local computer. The command output contains the URL for the dataset and an example ssh command to forward ports:

```
~/blobtoolkit/blobtools2/blobtools view --remote ~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME

    Initializing viewer |███████████████████████████████████████| 15/15 seconds
    Open dataset at http://localhost:8002/view/dataset/ASSEMBLY_NAME/blob?
    For remote access use:
        ssh -L 8002:127.0.0.1:8002 -L 8000:127.0.0.1:8000 username@remote_host
```

Host a dataset on your local computer and automatically open it in FireFox:

```
~/blobtoolkit/blobtools2/blobtools view --interactive ~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME
```

For more control of hosting parameters, use `blobtools host`:

```
~/blobtoolkit/blobtools2/blobtools host --port 8081 \
                                        --api-port 8001 \
                                        --hostname localhost \
                                        --viewer ~/blobtoolkit/viewer \
                                        ~/BTK_TUTORIAL/DATASETS

    Starting BlobToolKit API on port 8001 (pid: 76493)
    Starting BlobToolKit viewer on port 8081 (pid: 76495)
    Visit http://localhost:8081 to use the interactive BlobToolKit Viewer.
```

For even more control see the Viewer tutorial _Hosting a local instance_ for information on how to configure and run the Viewer directly.

#### `blobtools view`

- `--remote`
- `--interactive`
