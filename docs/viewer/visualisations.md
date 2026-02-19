---
layout: default
title: Visualisations
---

# Interpreting Visualisations

Guide to understanding BlobToolKit visualizations and plots.

## Blob Plot

The blob plot is the primary visualization in BlobToolKit:

- **X-axis**: GC content (percentage)
- **Y-axis**: Read coverage (depth)
- **Bubble size**: Contig length
- **Color**: Taxonomic assignment

### Reading the Plot

- **Dense clusters** indicate main organism groups
- **Outliers** may represent contamination or unusual sequences
- **Small bubbles** are short contigs (less reliable)

## Snail Plot

The snail plot shows cumulative contig coverage:

- **Spiral outward** from center represents ranked contigs by coverage
- **Color** indicates taxonomy
- **Useful for** identifying coverage distribution and anomalies

## Filtering by Taxonomy

Use the interactive legend to:

- Show/hide taxonomic groups
- Filter contigs by selection
- Export subsets

## Quality Metrics

Key metrics displayed:

- **N50/L50**: Contig length distribution
- **Total length**: Assembly size
- **GC content**: Overall genome GC%
- **Coverage**: Read depth statistics

## Common Patterns

- **Horizontal line at top**: Highly covered contigs (possible duplicates)
- **Scattered points**: Low coverage contigs (possible errors)
- **Multiple clusters**: Indicates mixed samples or contamination

## Next Steps

For command-line operations, see [BlobTools Commands](../command-line/commands).
