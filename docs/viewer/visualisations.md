---
layout: default
title: Visualisations
---

# Visualisations

Each BlobToolKit dataset is represented by several interactive views. All views
update automatically when filters or settings change.

## Changing Views

Available plot types are listed in the **Settings** menu (top of screen). Click
a view name to switch. Any view can be exported as SVG or PNG using the buttons
at the top right.

![Settings menu showing available views]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-09.02.4-1024x831.jpg' | relative_url }})

---

## Blob Plot

The default view when coverage data is present. Three shape variants are
available:

### Square-binned (default)

Contigs are binned along the x- and y-axes. Each square is scaled (square-root
scale) to represent the sum of scaffold lengths for each phylum in that bin.
Click any bin to select all scaffolds it contains — selected regions are
highlighted in pink in the plot and in the Summary panel.

- **X-axis**: mean GC content
- **Y-axis**: mean per-base coverage (log scale, broken at 0.01)
- **Colour**: taxonomic phylum (or chosen rank)
- **Size**: total span of scaffolds in bin

![Square-binned blob plot]({{ '/assets/img/viewer/ACVV01.blob_.square-1024x1024.png' | relative_url }})

Clicking bins highlights the selection in pink:

![Selected bins highlighted in pink]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-09.49.1-1024x831.jpg' | relative_url }})

Use a **log scale function** to restore prominence to small cobionts that a
square-root scale would de-emphasise:

![Log-scaled square-binned blob plot]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-09.43.36-1024x823.jpg' | relative_url }})

### Hex-binned

Same as square-binned but using hexagonal bins. Can have aesthetic advantages
when features are not axis-aligned.

![Hex-binned blob plot]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.04.08-1024x831.jpg' | relative_url }})

### Circle

One circle per scaffold — most similar to legacy BlobTools output. Drawing
individual shapes has performance implications for large datasets; see
[Performance](#performance) below.

![Circle blob plot]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.04.16-1024x831.jpg' | relative_url }})

---

## Cumulative View

Shows assembly length as scaffolds are added longest-to-shortest. The curve
shape reflects assembly contiguity. Per-phylum curves are overlaid.

![Cumulative view]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.04.43-1024x831.jpg' | relative_url }})

**Stacked variant**: per-phylum curves are normalised to share the same start
and end points as the overall curve — changes in direction highlight scaffold
length distribution across phyla.

![Cumulative stacked view]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.04.54-1024x831.jpg' | relative_url }})

Toggle the **curve origin scale** icon to rescale to the current _filtered_
dataset rather than the complete assembly.

---

## Snail Plot

Conveys multiple assembly metrics in a single graphic:

- Scaffolds arranged size-order clockwise from the outside of the central
  spiral, starting with the longest (red).
- Dark and light orange arcs mark **N50** and **N90**.
- Central light grey spiral shows cumulative scaffold count; white lines at
  each order of magnitude.
- Outer band: mean, max and min **GC vs AT** content at 0.1% intervals; white
  gaps reflect the proportion of Ns.
- **BUSCO scores** shown in the upper-right corner if available.

![Snail plot]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.05.32-1024x831.jpg' | relative_url }})

Hovering over any point shows the N*x* count and length in the lower-right
corner:

![Snail plot hover tooltip]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-11.23.33-1-1024x831.jpg' | relative_url }})

Click the scale boxes in the lower-left to set a fixed circumferential or
radial scale — useful for direct comparison of related assemblies:

![Snail plot with adjusted scale]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-11.24.32-1024x831.jpg' | relative_url }})

Toggle the **snail origin scale** icon to rescale to the filtered dataset.

---

## BUSCO View

BUSCO analyses are shown as tables and summary plots. For multi-lineage runs,
lineages are ordered most-specific first. Download buttons export all results
as JSON or CSV. Checkboxes in the left column select all scaffolds containing
BUSCOs in that subset.

![BUSCO view]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.05.07-1024x831.jpg' | relative_url }})

---

## Detail View

A formatted summary of dataset metadata with links to ENA, NCBI and Wikipedia
(for public/INSDC datasets). The full metadata can be downloaded as JSON. For
datasets produced with the BlobToolKit pipeline, the pipeline version is
included for reproducibility.

![Detail view]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.05.21-1024x831.jpg' | relative_url }})

---

## Table View

One row per scaffold with scores for all active fields. Click column headers to
sort. Export all data as CSV using the "csv" button. The Settings panel can be
hidden by clicking its tab to maximise table width.

![Table view]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-16.27.3-1024x831.jpg' | relative_url }})

Clicking the coloured square in the categories column for any scaffold opens a
**hit distribution plot** showing BLAST/Diamond results as a cumulative bitscore
plot coloured by taxon. Hover over individual hits for bitscore, position and
accession; click a hit to open its public database entry.

![Hit distribution plot]({{ '/assets/img/viewer/Screenshot-2019-07-29-at-16.33.15-1024x345.png' | relative_url }})

---

## Adjusting Plot Settings

Open the **Settings** menu to access plot-specific controls.

### Blob plot settings

| Setting                 | Description                                                                                                               |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| **Shape**               | Switch between square-binned, hex-binned and circle variants.                                                             |
| **Resolution**          | Number of bins per axis. Higher resolution separates overlapping blobs; lower makes region selection easier.              |
| **Reducer function**    | How scaffold lengths in each bin are combined: sum (default), min, max, count, mean.                                      |
| **Scale function**      | How values are scaled: square-root (default, area ∝ value), log (emphasises small values), linear.                        |
| **Scale factor**        | Maximum bin display size relative to bin width. Increase to make blobs more prominent and allow adjacent bins to overlap. |
| **Palettes**            | Choose a colour palette. Select **custom** then click swatches to assign per-category colours.                            |
| **PNG resolution (px)** | Width in pixels of exported PNG — useful for publication-quality output.                                                  |
| **Static threshold**    | Maximum dataset size for interactive display. Datasets above this size show pre-rendered images.                          |

Setting the **reducer function** to _count_ and **scale function** to _log_
highlights scaffold count per bin rather than span:

![Count + log scale]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-15.26.21-1024x831.jpg' | relative_url }})

Increasing the **scale factor** makes blob regions more prominent:

![Increased scale factor]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-15.47.54-1024x831.jpg' | relative_url }})

Changing the **resolution** adjusts how many bins the axes are divided into.
Lower values make selection easier; higher values reveal fine-scale patterns:

|                             Resolution 40                              |  Resolution 30   |                             Resolution 20                              |
| :--------------------------------------------------------------------: | :--------------: | :--------------------------------------------------------------------: | ---------------- | ---------------------------------------------------------------------- | ---------------- |
| ![res 40]({{ '/assets/img/viewer/ACVV01.blob\_.square-6-1024x1024.png' | relative_url }}) | ![res 30]({{ '/assets/img/viewer/ACVV01.blob\_.square-5-1024x1024.png' | relative_url }}) | ![res 20]({{ '/assets/img/viewer/ACVV01.blob\_.square-4-1024x1024.png' | relative_url }}) |

Selecting the _custom_ **palette** and clicking swatches lets you assign
publication-ready category colours:

![Custom palette editor]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-15.58.31-1024x823.jpg' | relative_url }})

### Snail / cumulative settings

The **snail origin scale** and **curve origin scale** icons rescale axes to
the current filtered dataset rather than the complete assembly — helpful when a
large portion of the assembly has been filtered out:

![Snail filtered (unscaled)]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-16.27.35-1024x823.jpg' | relative_url }})

![Snail filtered (origin scaled)]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-16.27.45-1024x831.jpg' | relative_url }})

---

## Performance

For datasets with more than 100,000 scaffolds the viewer shows pre-rendered
static images by default. To view large datasets interactively:

1. Increase the **static threshold** setting in the Settings menu, **or**
2. Click the **interactive** button at the top of the Settings menu.

Additional tips for large datasets:

- Use the square-binned or hex-binned blob plot rather than circles.
- Reduce the **resolution** to decrease the number of rendered bins.
- Reduce the **scale factor** to minimise rendering overlap.

### Plot graphics threshold

The **plot graphics threshold** (default 10,000) controls when circle plots
switch from SVG to HTML5 canvas rendering. Above the threshold the browser
draws the plot as a bitmap, which uses less memory but cannot be exported as an
image file. To force SVG rendering, increase the threshold above the number of
scaffolds in the dataset, or click the SVG icon in the Settings menu:

![Plot graphics threshold setting]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.27.36-1024x831.jpg' | relative_url }})

### Circle limit

The **circle limit** (default 100,000) reduces render time by omitting
no-hit circles until the total number of circles would exceed the limit. The
threshold can be raised to show all circles or lowered to keep the plot
responsive:

![Circle limit setting]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.28.19-1024x831.jpg' | relative_url }})

### Static threshold

On the public viewer, datasets above the **static threshold** (default
100,000 scaffolds) are presented as pre-rendered static images. A warning is
shown when static rendering is active; increase the threshold or click the
interactive button to enable full interactivity:

![Static threshold warning]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-10.28.47-1024x831.jpg' | relative_url }})

---

## Changing Plot Axes

Blob plots show GC vs. coverage by default. Any available field can be
assigned to either axis using the _Filters_ menu. Active fields are shown with
a pink header and a preview histogram. Icons in the header indicate the current
axis assignment (icons with a dark background are active). The default GC /
coverage / length (x/y/z) configuration looks like this:

![Default filters menu with axis icons]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-12.09.28-1024x831.jpg' | relative_url }})

#### Coverage–coverage plots

Datasets with more than one coverage library can be plotted as cov–cov plots to
highlight scaffolds with divergent coverage ratios (potential contaminants).
Activate the second coverage field by clicking its header:

![Second coverage field activated]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-12.11.24-1024x831.jpg' | relative_url }})

Then click the x-axis icon in that field's header. Most scaffolds fall on the
diagonal; those far from it have very different relative coverage in the two
libraries:

![Coverage–coverage blob plot]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-12.11.38-1024x831.jpg' | relative_url }})

#### Changing category rank

Category colours are set by the `bestsumorder_phylum` field by default. Scroll
down the Filters menu to see the category previews:

![Category field at phylum rank]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-13.41.07-1024x823.jpg' | relative_url }})

To change the rank to family, activate the `bestsumorder_family` field and click
its category icon:

![Category field switched to family rank]({{ '/assets/img/viewer/Screenshot-2019-07-30-at-13.45.57-1024x831.jpg' | relative_url }})

---

## Filtering Assemblies

Assemblies can be filtered interactively by setting value ranges, toggling
taxonomic categories, or applying a selection-based filter — all from the
_Filters_ menu.

### Slider-based filtering

Drag either slider on a preview histogram to set an inclusive range. The plot
updates when you release the slider. The examples below use log-scaled counts
to make the effect more visible (see [Adjusting plot settings](#settings)):

![Opening filters and dragging slider]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.20.56-1024x831.jpg' | relative_url }})

The current view updates to reflect the filter — here, high- and low-GC
scaffolds have been excluded:

![GC filter applied]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.21.39-1024x831.jpg' | relative_url }})

Filters carry over when switching between views:

|                          Cumulative                           |      Snail       |                      Blob (circle)                       |
| :-----------------------------------------------------------: | :--------------: | :------------------------------------------------------: | ---------------- | ---------------------------------------------------------------- | ---------------- |
| ![]({{ '/assets/img/viewer/ACVV01.cumulative-3-1024x1024.png' | relative_url }}) | ![]({{ '/assets/img/viewer/ACVV01.snail-2-1024x1024.png' | relative_url }}) | ![]({{ '/assets/img/viewer/ACVV01.blob\_.circle-2-1024x1024.png' | relative_url }}) |

Type exact values into the min/max boxes to set precise boundaries, for example
setting a minimum length of 5,000 bp:

![Minimum length set to 5000]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.22.31-1024x831.jpg' | relative_url }})

Click the **invert** icon in any field header to exclude rather than include
scaffolds matching the filter:

![Inverted length filter]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.23.12-1024x831.jpg' | relative_url }})

### Category-based filtering

Click the coloured swatches at the top of a category preview histogram to toggle
individual taxa on or off. The example below excludes "no-hit" sequences:

![Category filter — no-hit excluded]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.24.37-1024x831.jpg' | relative_url }})

### Selection-based filtering

First create a selection on the blob plot (see [Using Selections](#using-selections)),
then click the **selection** filter header in the Filters menu to restrict the
view to that selection:

![Arbitrary selection drawn on blob plot]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.50.40-1024x831.jpg' | relative_url }})

![Selection filter activated]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.51.30-1024x831.jpg' | relative_url }})

Selection filters can also be inverted and combined with other filters:

![Inverted selection combined with coverage filter]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.51.39-1024x831.jpg' | relative_url }})

![Combined filter result]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-09.53.28-1024x831.jpg' | relative_url }})

---

## Using Selections

Selections let you highlight specific scaffolds and track them across different
views. Selections can be made in blob, BUSCO or table view.

### Making a BUSCO-based selection

Open the BUSCO view from the Settings menu. Click the checkbox to the left of
**Complete** in the BUSCO results table to select all scaffolds with a complete
BUSCO gene. Other checkboxes showing a half-filled state indicate that only some
scaffolds in those categories are included:

![BUSCO view — opening the settings]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.20.52-1024x823.jpg' | relative_url }})

![Complete BUSCO genes selected]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.25.09-1024x831.jpg' | relative_url }})

Open the **Summary** tab to see statistics for the selected scaffolds (shown in
pink). In this example the selection accounts for ~35% of total assembly span,
with over 60% assigned to Arthropoda, but only 3.7% of scaffold count:

![Summary for selected scaffolds]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.46.44-1024x831.jpg' | relative_url }})

### Tracking selections across views

Switch to blob view in the Settings menu. Bins containing any selected scaffolds
are outlined in pink. Most cluster around the Arthropoda peak, but four bins at
higher GC form a separate cluster:

![Selected scaffolds visible on blob plot]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.34.38-1024x831.jpg' | relative_url }})

### Investigating outliers with filters

Set a minimum GC filter of 0.55 in the Filters menu to isolate the high-GC
selected scaffolds:

![GC filter isolating high-GC selected bins]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.35.03-1024x831.jpg' | relative_url }})

Switch to table view to inspect the individual scaffolds:

![Table view with filtered selection]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.35.27-1024x831.jpg' | relative_url }})

Click the checkbox column header to sort selected scaffolds to the top of the
list:

![Selected scaffolds sorted to top]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.36.50-1024x831.jpg' | relative_url }})

All selected scaffolds are labelled Proteobacteria — likely a _Wolbachia_
endosymbiont. Click a coloured square in the Categories column to see the
underlying BLAST hit distribution and links to NCBI gene records:

![Hit distribution for a selected scaffold]({{ '/assets/img/viewer/Screenshot-2019-07-31-at-10.39.27-1024x831.jpg' | relative_url }})

---

## Next Steps

- [Local Hosting](local-hosting) — start a local viewer instance.
- [BlobTools Commands](../command-line/commands) — filtering and exporting
  datasets on the command line.
