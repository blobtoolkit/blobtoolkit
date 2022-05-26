#!/usr/bin/env python3
# pylint: disable=too-many-locals,too-many-branches,too-many-arguments,too-many-statements,too-many-nested-blocks

"""Parse text files into fields."""

import math
import pathlib
import re
from collections import defaultdict

from ..lib import file_io
from .field import Array
from .field import Category
from .field import MultiArray
from .field import Variable


def parse_header_row(delimiter, header_row, columns):
    """Parse header row."""
    head = re.split(delimiter, header_row)
    for i, col in enumerate(columns):
        try:
            index, name = col.split("=")
            if not index.isdigit():
                try:
                    index = head.index(index)
                    columns[i] = "%d=%s" % (index + 1, name)
                except ValueError:
                    pass
        except ValueError:
            if col.isdigit():
                try:
                    columns[i] = "%s=%s" % (col, head[int(col) - 1])
                except IndexError:
                    exit("ERROR: column index out of range '%s'" % col)
            else:
                try:
                    index = head.index(col)
                    columns[i] = "%s=%s" % (str(index + 1), col)
                except ValueError:
                    pass
    if not columns:
        return head
    return columns


def map_fields(delimiter, first_row, columns):
    """Set column headers and field types."""
    sample_data = re.split(delimiter, first_row)
    width = len(sample_data)
    cols = []
    headers = {}
    types = {}
    for i, col in enumerate(columns):
        if col:
            try:
                index, name = col.split("=")
                if index.isdigit() and width >= int(index):
                    cols.append(int(index) - 1)
                    headers[int(index) - 1] = name
                else:
                    exit("ERROR: column index out of range '%s'" % col)
            except ValueError:
                if col.isdigit():
                    exit("ERROR: no name specified for column %s" % col)
                else:
                    cols.append(i)
                    headers[i] = col
    for col in cols:
        if re.search("identifier", headers[col], re.IGNORECASE):
            types[headers[col]] = "Identifier"
        else:
            datum = sample_data[col]
            try:
                float(datum)
                types[headers[col]] = "Variable"
            except ValueError:
                types[headers[col]] = "Category"
    return cols, headers, types, width


def parse_rows(delimiter, lines, width, no_array, cols, types, headers):
    """Parse rows and test for duplicate identifiers."""
    ids = set()
    rows = []
    id_rows = defaultdict(list)
    index = 0
    array = False
    for line in lines:
        line = line.replace('"', "")
        data = re.split(delimiter, line)
        if len(data) < width:
            continue
        row = {}
        for col in cols:
            if types[headers[col]] == "Identifier":
                if data[col] in ids:
                    if no_array:
                        exit(
                            "ERROR: found multiple instances of Identifier '%s'"
                            % data[col]
                        )
                    array = True
                ids.add(data[col])
                id_rows[data[col]].append(index)
                index += 1
            else:
                row[headers[col]] = data[col]
        rows.append(row)
    return rows, id_rows, array, ids


def rows_to_results(rows, id_rows, types, array, field_name):
    """Make fields from rows."""
    field_names = [type for type in types.keys() if types[type] != "Identifier"]
    if field_name:
        results = {field_name: {}}
    else:
        results = {name: {} for name in field_names}
    for ident, indices in id_rows.items():
        if field_name:
            if len(field_names) > 1:
                results[field_name][ident] = []
        elif array:
            for name in field_names:
                results[name][ident] = []
        for index in indices:
            row = rows[index]
            if field_name:
                data = [row[name] for name in field_names]
                if len(field_names) == 1:
                    data = data[0]
                if array:
                    results[field_name][ident].append(data)
                else:
                    results[field_name][ident] = data
            else:
                for name in field_names:
                    if array:
                        results[name][ident].append([row[name]])
                    else:
                        results[name][ident] = row[name]
    return results


def results_to_fields(results, types, cols, headers, text_file, delimiter, identifiers):
    """Convert results to fields."""
    fields = []
    field_types = {
        "Variable": Variable,
        "Category": Category,
        "Array": Array,
        "MultiArray": MultiArray,
    }
    array_type = "array"
    for col in cols:
        if types[headers[col]] == "Identifier":
            id_column = col
        else:
            if array_type and types[headers[col]] != array_type:
                array_type = "mixed"
                break
            array_type = "string" if types[headers[col]] == "Category" else "float"
    for key, values in results.items():
        ident, sample = next(iter(values.items()))
        blank = "NA" if types[headers[col]] == "Category" else 0
        if not identifiers.validate_list(list(results[key].keys())):
            print(
                "WARN: Contig names in the text file did not match dataset identifiers."
            )
        kwargs = {
            "meta": {
                "field_id": key,
                "name": key,
                "preload": False,
                "active": False,
                "file": text_file,
                "id_column": id_column,
                "delimiter": delimiter,
                "datatype": array_type,
            },
            "parents": [],
        }
        if isinstance(sample, list):
            blank = []
            array_headers = []
            index = -1
            if key in types:
                array_headers.append(key)
                if types[key] == "Category" and "category_slot" not in kwargs:
                    kwargs.update({"category_slot": 0})
            else:
                for col in cols:
                    if types[headers[col]] != "Identifier":
                        index += 1
                        array_headers.append(headers[col])
                        if (
                            types[headers[col]] == "Category"
                            and "category_slot" not in kwargs
                        ):
                            kwargs.update({"category_slot": index})
            kwargs.update({"headers": array_headers})
            if isinstance(sample[0], list):
                field_type = field_types["MultiArray"]
                kwargs.update({"type": "multiarray"})
            else:
                field_type = field_types["Array"]
                kwargs.update({"type": "array"})
        else:
            field_type = field_types[types[key]]
            kwargs.update({"type": types[key].lower()})

        if kwargs["type"] == "variable":
            vals = []
            is_float = False
            for ident in identifiers.values:
                value = results[key][ident] if ident in results[key] else blank
                vals.append(value)
                if kwargs["type"] == "variable" and not is_float:
                    try:
                        int(value)
                    except ValueError:
                        is_float = True
            min_max = [math.inf, -math.inf]
            values = []
            for value in vals:
                value = float(value) if is_float else int(value)
                values.append(value)
                min_max = [min(min_max[0], value), max(min_max[1], value)]
            kwargs["meta"].update({"range": min_max})
            if is_float:
                kwargs["meta"].update({"datatype": "float"})
            else:
                kwargs["meta"].update({"datatype": "integer"})
            if min_max[0] < 0 or min_max[0] > min_max[1] / 1000:
                kwargs["meta"].update({"scale": "scaleLinear"})
            else:
                kwargs["meta"].update({"scale": "scaleLog"})
                if min_max[0] == 0:
                    kwargs["meta"].update({"clamp": 0.01})
        else:
            values = [
                results[key][ident] if ident in results[key] else blank
                for ident in identifiers.values
            ]
            if kwargs["type"] == "category":
                kwargs["meta"].update({"datatype": "string"})
        field = field_type(key, values=values, **kwargs)
        fields.append(field)
    return fields


def set_delimiter(delimiter, *, sample=None):
    """Set text delimiter."""
    if delimiter == "whitespace":
        if sample is not None:
            if "\t" in sample:
                return re.compile(r"\t")
        return re.compile(r"\s+")
    else:
        return re.compile(r"%s" % delimiter)


def parse_text(text_file, delimiter, columns, header, no_array, identifiers):
    """Parse text file into Category and/or Variable fields."""
    try:
        text_file, field_name = text_file.split("=")
    except ValueError:
        field_name = False
    data = file_io.read_file(text_file)
    lines = data.split("\n")
    delimit = set_delimiter(delimiter, sample=lines[0])
    if columns:
        columns = columns.split(",")
    else:
        columns = []
    if header:
        header_row = lines.pop(0).replace('"', "")
        columns = parse_header_row(delimit, header_row, columns)
    cols, headers, types, width = map_fields(
        delimit, lines[0].replace('"', ""), columns
    )
    rows, id_rows, array = parse_rows(
        delimit, lines, width, no_array, cols, types, headers
    )[:3]
    # if not identifiers.validate_list(list(ids)):
    #     exit('ERROR: contig names in the text file did not match dataset identifiers.')
    results = rows_to_results(rows, id_rows, types, array, field_name)
    fields = results_to_fields(
        results, types, cols, headers, text_file, delimiter, identifiers
    )
    # meta = {'file': text_file}
    return fields
    # results = defaultdict(list)
    # for line in lines:
    #     if header:
    #         row = re.split(' +', line)
    #         if len(row) > 1:
    #             if row[1].startswith('v.'):
    #                 meta.update({'version': row[1]})
    #             elif row[1] == 'Mode:':
    #                 meta.update({'mode': row[2]})
    #                 meta.update({'field_id': "trnascan_%s" % row[2].lower()})
    #             elif row[1].startswith('------'):
    #                 header = False
    #     else:
    #         row = re.split(r' +|\t', line)
    #         if len(row) == 9:
    #             results[row[0]].append([row[4], row[5]])
    # if not identifiers.validate_list(list(results.keys())):
    #     raise UserWarning('Contig names in the tRNAScan file did not match dataset identifiers.')
    # values = [results[id] if id in results else [] for id in identifiers.values]
    # trnascan_field = MultiArray(meta['field_id'],
    #                             values=values,
    #                             meta=meta,
    #                             headers=('tRNA_type', 'Anticodon'),
    #                             parents=['children']
    #                             )
    # return trnascan_field


def apply_filter(ids, text_file, **kwargs):
    """Filter Text file."""
    suffix = kwargs["--suffix"]
    path = pathlib.Path(text_file)
    outfile = str(path.parent / (path.stem + "." + suffix + path.suffix))
    data = file_io.read_file(text_file)
    lines = data.split("\n")
    delimiter = kwargs["--text-delimiter"]
    delimit = set_delimiter(delimiter, sample=lines[0])
    id_col = int(kwargs["--text-id-column"]) - 1
    output = []
    if kwargs["--text-header"]:
        header_row = lines.pop(0)
        header_row.rstrip()
        output.append(header_row)
    for line in lines:
        line = line
        row = re.split(delimit, line.replace('"', ""))
        try:
            if row[id_col] in ids:
                output.append(line)
        except IndexError:
            output.append(line)
    file_io.write_file(outfile, output, plain=True)


def parse(files, **kwargs):
    """Parse all text files."""
    parsed = []
    for file in files:
        try:
            fileData = parse_text(
                file,
                delimiter=kwargs["--text-delimiter"],
                columns=kwargs["--text-cols"],
                header=kwargs["--text-header"],
                no_array=kwargs["--text-no-array"],
                identifiers=kwargs["dependencies"]["identifiers"],
            )
        except AttributeError:
            continue
        parsed = parsed + fileData
    return parsed


def parent():
    """Set standard metadata for text."""
    return []
