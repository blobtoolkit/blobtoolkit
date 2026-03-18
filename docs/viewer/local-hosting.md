---
layout: default
title: Local Hosting
---

# Local Hosting

To explore locally generated BlobDir datasets interactively you can start a
local viewer instance using the `blobtools host` command.

## Quick Start

```bash
blobtools host /path/to/datasets/
```

This starts both the API and the viewer client on their default ports. Open
`http://localhost:8080/view/all` in a browser to see your datasets.

## Port Configuration

```bash
blobtools host --port 8080 --api-port 8000 /path/to/datasets/
```

The `--port` flag sets the client port (default 8080). The `--api-port` flag
sets the API port (default 8000).

## Connecting from a Remote Server

If BlobTools is running on a remote server over SSH, forward both ports before
opening your browser:

```bash
ssh -L 8000:127.0.0.1:8000 -L 8080:127.0.0.1:8080 user@remote-host
```

Then visit `http://localhost:8080/view/all` locally.

## Environment Variable Configuration

For persistent settings, copy the provided `.env.dist` to `.env` in the viewer
directory and edit it:

```bash
cp .env.dist .env
```

Key variables:

| Variable                | Default                | Description                                                  |
| ----------------------- | ---------------------- | ------------------------------------------------------------ |
| `BTK_FILE_PATH`         | —                      | **Required.** Path to directory containing BlobDir datasets. |
| `BTK_CLIENT_PORT`       | `8080`                 | Port for the viewer client.                                  |
| `BTK_API_PORT`          | `8000`                 | Port for the API.                                            |
| `BTK_HOST`              | `localhost`            | Hostname if accessing via a name other than localhost.       |
| `BTK_ORIGINS`           | _(localhost variants)_ | CORS allowed origins for API connections.                    |
| `BTK_DATASET_TABLE`     | `false`                | Display datasets in a table on the main landing page.        |
| `BTK_USE_DEFAULT_LINKS` | `true`                 | Include links to external resources (ENA, NCBI, Wikipedia).  |
| `BTK_CIRCLE_LIMIT`      | `100000`               | Maximum records for which individual circles are plotted.    |
| `BTK_NOHIT_THRESHOLD`   | `1000000`              | Maximum records for which no-hit data is plotted.            |

## Multiple Datasets

Point `BTK_FILE_PATH` (or the command-line path argument) at a **directory
containing multiple BlobDirs** rather than a single BlobDir. All valid BlobDirs
found will be listed on the viewer landing page.

```
datasets/
  GCA_000001405/
  GCA_000950515/
  GCA_949316315/
```

```bash
blobtools host datasets/
```

## Performance Notes

- For datasets with more than 100,000 scaffolds the viewer will display
  pre-rendered static images by default. To enable interactive viewing for
  large datasets, increase BTK_CIRCLE_LIMIT or click the "interactive" button
  in the Settings menu.
- Reducing the blob plot **resolution** setting (number of bins on each axis)
  can improve responsiveness for large datasets.

## Public vs. Local Instances

We host a public instance of the viewer at
[blobtoolkit.genomehubs.org/view](https://blobtoolkit.genomehubs.org/view)
where thousands of publicly available datasets can be explored. Both the public
and local instances share all core features.

### Public instance home page

The public home page includes taxonomic progress bars tracking how many
INSDC-registered assemblies have been analysed with the BlobToolKit pipeline:

![Public viewer home page with taxonomy progress bars]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-12.04.1-1024x831.jpg' | relative_url }})

Search results on the public instance are shown as a rich summary-statistics
table (pre-computed via the CLI). The table can be sorted, filtered, and
customised via the **Customise table** link, and exported as a CSV file:

![Public instance search results table]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-12.08.51-1024x831.jpg' | relative_url }})

### Local instance home page

The default local home page shows a simple search bar on the left. This layout
can be altered to resemble the public instance by setting the `BTK_*`
environment variables described above:

![Local viewer home page]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-12.05.22-1024x831.jpg' | relative_url }})

Local search results are displayed as a plain list showing the dataset ID, GCA
accession, taxon name and number of records (pre-computed summary statistics are
not shown unless the environment variables are configured):

![Local instance search results list]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-14.39.5-1024x831.jpg' | relative_url }})

---

## Searching Available Datasets

### Using the search box

Type any search term into the search box at the top of the page.
Autocomplete suggestions appear as you type based on dataset metadata. Here,
searching for "Nematoda" shows suggestions after typing "Ne":

![Search box with autocomplete suggestions]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-15.14.0-1024x831.jpg' | relative_url }})

Select a suggestion or press Enter to see the matching datasets. The rich table
view below is available on the public instance where pre-computed summary
statistics exist:

![Search results in summary-statistics table]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-15.14.50-1024x831.jpg' | relative_url }})

Click any column header to sort the results. For example, clicking "sequences"
orders by scaffold count:

![Results sorted by sequence count]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-15.15.16-1024x831.jpg' | relative_url }})

Click **Customise table** to show or hide columns. Use the **csv** link below
the table to download the results:

![Customise table panel]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-15.35.0-1024x831.jpg' | relative_url }})

Scroll below the table to see the taxonomy browser automatically expanded to the
rank that best matches your search term (the matching taxon is marked with a
vertical pink bar):

![Taxonomy browser expanded to matching rank]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-15.42.51-1024x831.jpg' | relative_url }})

### Using the taxonomy browser

The taxonomy browser lets you navigate INSDC assemblies at eight common ranks
(superkingdom → species). Progress bars show how many datasets are available on
the public instance for each taxon. Click parenthesised numbers to expand or
collapse nodes:

![Taxonomy browser navigation]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-15.55.5-1024x831.jpg' | relative_url }})

Click a taxon name to search for all matching datasets:

![Taxon name search results]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-15.56.11-1024x831.jpg' | relative_url }})

---

## Next Steps

- [Visualisations Guide](visualisations) — interpreting blob, snail, cumulative
  and table views.
- [Creating a dataset](../command-line/tutorials) — generating a BlobDir to
  host.
