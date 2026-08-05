---
title: Adding text files
date: 2020-02-12
---

Files can be imported from generic text files, with the ability to specify column separators and map columns to field names. This provides flexibility to view and filter a wide range of additional analysis outputs beyond those with dedicated parsers.

### Example

Text file import can be from any text file. This example just uses GC proportion and length of the 10 longest contigs in the _Strongyloides venezuelensis_ assembly

cd ~/BTK_TUTORIAL/FILES

nano ASSEMBLY_NAME.values.csv

"id","gc","length"
"LM524968.1",0.2553,5851329
"LM524969.1",0.2523,3462877
"LM524970.1",0.2523,2731838
"LM524971.1",0.2513,2360143
"LM524972.1",0.2711,1484661
"LM524973.1",0.2619,1243793
"LM524974.1",0.2591,1114134
"LM524975.1",0.2558,1106993
"LM524976.1",0.2374,1082826
"LM524977.1",0.2577,1069219

This file can be imported to add two new fields to the dataset:

```
~/blobtoolkit/blobtools2/blobtools add \
    --text ASSEMBLY_NAME.values.csv \
    --text-delimiter ',' \
    --text-cols id=identifiers,gc=gc_proportion,length=contig_length \
    --text-header \
    ~/BTK_TUTORIAL/DATASETS/ASSEMBLY_NAME
```

### Configuration

There are several configuration options when adding text files using `blobtools add`:

- `[--text](#cov)` – Text file to import. _\[_**_Required_**_\]_
- `[--delimiter](#delimiter)` – Text file delimiter. _\[Default: whitespace\]_
- `[--text-cols](#text-cols)` – Comma separated list of _index\[=field_name\]_ or _header\[=field_name\]_. _\[_**_Required_**_\]_
- `[--text-header](#text-header)` – Flag to indicate if first row of text file contains field names. _\[Default: false\]_
- `[--text-no-array](#text-no-array)` – Flag to prevent fields in files with duplicate identifiers being loaded as array fields. _\[Default: false\]_
- `[--replace](#replace)` – Replace existing fields if present. _\[Default: false\]_

#### `--text`

A text file can be specified using the `--text` flag. As with other parsers, multiple files may be specified, however, it would be necessary for field names to be specified using unique headers in the different files due to the mapping of field names to columns. Alternatively, it is possible to specify `=FIELDNAME` after the filename to load all specified columns into a single _array_ or _multiarray_ type field, with multiple values (or multiple sets of values) per identifier.

**Usage:**

```
blobtools add \
    --text ASSEMBLY_NAME.results.txt \
    --text ASSEMBLY_NAME.results2.txt=results2 \
    ...
```

#### `--delimiter`

The text file column separator can be specified using the `--delimiter` flag. The default value is `whitespace`, which spilts rows on any whitespace character (e.g. SPACE or TAB). Alternative delimiters can use [Python regular expression syntax](https://docs.python.org/3/howto/regex.html), e.g. use `--delimiter '\t'` for tab delimited files.

**Usage:**

```
blobtools add \
    ...
    --delimiter '\t' \
    ...
```

#### `--text-cols`

For any text file, it is necessary to specify which column contains the contig identifiers (which must match those already imported into the dataset from sequence headers) and at least one column containing values to be assigned to a dataset field. Column specification can be based on column indices or (if `[--text-header](#text-header)` is set) column names. Field names can be set by adding `=FIELD_NAME` after the column index/header (field names must be set this way if `[--text-header](#text-header)` is not set). The datatype for each field (category or variable) will be detected automatically during import. If the file contains duplicate identifiers and `[--text-no-array](#text-no-array)` is not set, _array_ fields will be created with multiple values per identifier.

**Usage:**

```
blobtools add \
    --text-cols 1=identifiers,2,3=score,total=total_score \
    ...
```

#### `--text-header`

Set the `--text-header` flag if a text file contains a header row. if present, the header row can be used to determine field names.

**Usage:**

```
blobtools add \
    --text-header \
    ...
```

#### `--text-no-array`

Set the `--text-no-array` flag to prevent files with duplicated identifiers being loaded as _array_ fields. Typically used when a text file is not expected to contain duplicate identifiers.

**Usage:**

```
blobtools add \
    --text-no-array \
    ...
```

#### `--replace`

If a `blobtools add` command would overwrite an existing field, the default behaviour is to issue a warning and not replace the existing field. To change this behaviour and allow existing fields to be overwritten, set the `--replace` flag.

**Usage:**

```
blobtools add \
    ...
    --replace \
    ...
```
