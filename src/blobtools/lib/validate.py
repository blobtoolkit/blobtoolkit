#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-statements, too-many-locals, W0603, W0703

"""
Validate a BlobDir.

Usage:
    blobtools validate [--basic] [--example] DIRECTORY

Arguments:
    DIRECTORY             BlobDir directory.

Options:
    --basic    Only require basic metadata.
    --example  Validate example dataset.
"""

import copy
import os
import os.path
import re
import sys
import urllib

import fastjsonschema
import ujson
from docopt import docopt

from .version import __version__

PIDS = []


def load_json_file(file, schema=False):
    """Load a JSON file."""
    result = []
    error_offset = 7 if schema else 0
    try:
        with open(file, "r") as fh:
            data = ujson.load(fh)
            if "properties" in data:
                for property, value in data["properties"].items():
                    if "$ref" in value:
                        ref_file = os.path.join(os.path.dirname(file), value["$ref"])
                        if os.path.exists(ref_file):
                            value["$ref"] = "file://%s" % ref_file
            return data
    except IOError as err:
        show_error(err.strerror, 1 + error_offset)
    except ValueError as err:
        show_error(err, 2 + error_offset)


def url_scheme(url, path):
    """Treat local URLs as 'file://'."""
    if not urllib.parse(url).scheme:
        url = "file://" + os.path.join(path, url)
    return url


def show_error(error, step=0):
    """List errors then quit."""
    stages = [
        "during undefined stage",
        "while trying to read file",
        "while trying to parse JSON",
        "while testing raw metadata against schema",
        "while testing processed metadata against schema",
        "while generating data schema",
        "while testing data against generated schema",
        "while checking standard fields",
        "while trying to read schema file",
        "while trying to parse JSON schema",
    ]
    if error:
        stage = stages[step]
        print("NOT VALID: Error %s" % stage)
        quit(error)


def flatten_fields(fields, parent={}):
    """Flatten hierarchical field list."""
    field_list = []
    for field in fields:
        for key, value in parent.items():
            if key not in field and key not in ["children", "data"]:
                field[key] = value
        if "children" not in field:
            try:
                check_expected_field_properties(field)
            except AssertionError as err:
                show_error(err, 7)
            field_list.append(field)
        else:
            field_list += flatten_fields(field["children"], field)
        if "data" in field:
            field_list += flatten_fields(field["data"], field)
    return field_list


def check_expected_field_properties(field):
    """Test properties in expected fields."""
    field_id = field["id"]
    if field_id in ["gc", "length", "ncount"] or re.match(r".+_cov$", field_id):
        assert field["type"] == "variable", (
            "expected %s.type to be 'variable'" % field_id
        )
        assert field["range"][0] >= 0, "expected %s.range[0] to be >= 1" % field_id
        if field_id in ["length", "ncount"] or re.match(r".+read_cov$", field_id):
            assert field["datatype"] == "integer", (
                "expected %s.datatype to be 'integer'" % field_id
            )
        else:
            assert field["datatype"] == "float", (
                "expected %s.datatype to be 'float'" % field_id
            )
    elif re.match(r".+_busco", field_id):
        assert field["type"] == "multiarray", (
            "expected %s.type to be 'multiarray'" % field_id
        )
        assert field["datatype"] == "mixed", (
            "expected %s.datatype to be 'mixed'" % field_id
        )
        assert "Busco id" in field["headers"], (
            "expected to find 'Busco id' in %s.headers" % field_id
        )
        assert "Status" in field["headers"], (
            "expected to find 'Status' in %s.headers" % field_id
        )
        assert "category_slot" in field, (
            "expected to find 'category_slot' in %s" % field_id
        )
        message = (
            "expected %s.category_slot to match index of 'Status' in %s.headers"
            % (field_id, field_id)
        )
        assert field["category_slot"] == field["headers"].index("Status"), message
    elif re.match(r"bestsum.+", field_id):
        if re.match(r".+_cindex", field_id):
            assert field["type"] == "variable", (
                "expected %s.type to be 'variable'" % field_id
            )
            assert field["datatype"] == "integer", (
                "expected %s.datatype to be 'integer'" % field_id
            )
            assert field["range"][0] >= 0, "expected %s.range[0] to be >= 1" % field_id
        elif re.match(r".+_positions", field_id):
            assert field["type"] == "multiarray", (
                "expected %s.type to be 'multiarray'" % field_id
            )
            if re.match(r".+_.+_positions", field_id):
                assert field["datatype"] == "string", (
                    "expected %s.datatype to be 'string'" % field_id
                )
                assert "name" in field["headers"], (
                    "expected to find 'name' in %s.headers" % field_id
                )
                assert "category_slot" in field, (
                    "expected to find 'category_slot' in %s" % field_id
                )
                message = (
                    "expected %s.category_slot to match index of 'name' in %s.headers"
                    % (field_id, field_id)
                )
                assert field["category_slot"] == field["headers"].index("name"), message
                assert "linked_field" in field, (
                    "expected to find 'linked_field' in %s" % field_id
                )
            else:
                assert field["datatype"] == "mixed", (
                    "expected %s.datatype to be 'mixed'" % field_id
                )
                assert "taxid" in field["headers"], (
                    "expected to find 'taxid' in %s.headers" % field_id
                )
                assert "score" in field["headers"], (
                    "expected to find 'score' in %s.headers" % field_id
                )
                assert "start" in field["headers"], (
                    "expected to find 'start' in %s.headers" % field_id
                )
                assert "end" in field["headers"], (
                    "expected to find 'end' in %s.headers" % field_id
                )
                assert "subject" in field["headers"], (
                    "expected to find 'subject' in %s.headers" % field_id
                )
                assert "index" in field["headers"], (
                    "expected to find 'index' in %s.headers" % field_id
                )
        elif re.match(r".+_score", field_id):
            assert field["type"] == "variable", (
                "expected %s.type to be 'variable'" % field_id
            )
            assert field["datatype"] == "float", (
                "expected %s.datatype to be 'float'" % field_id
            )
            assert field["range"][0] >= 0, "expected %s.range[0] to be >= 1" % field_id
        else:
            assert field["type"] == "category", (
                "expected %s.type to be 'category'" % field_id
            )


def generate_data_schema(field, meta, data_schemas, data):
    """Modify default schema to match field properties."""
    data_schema = copy.deepcopy(data_schemas[field["type"]])
    schema_values = data_schema["properties"]["values"]
    schema_values["minItems"] = meta["records"]
    schema_values["maxItems"] = meta["records"]
    if field["type"] == "variable":
        schema_values["items"]["minimum"] = field["range"][0]
        schema_values["items"]["maximum"] = field["range"][1]
        value_type = "integer" if field["datatype"] == "integer" else "number"
        schema_values["items"]["type"] = value_type
    elif field["type"] == "category":
        try:
            schema_values["items"]["minimum"] = 0
            schema_values["items"]["maximum"] = len(data["keys"]) - 1
        except KeyError as err:
            show_error(err, 5)
    elif field["type"] == "multiarray":
        schema_values["items"]["items"]["minItems"] = len(field["headers"])
        schema_values["items"]["items"]["maxItems"] = len(field["headers"])
        if "category_slot" in field:
            try:
                schema_values["items"]["items"]["items"] = [
                    {} for x in field["headers"]
                ]
                category = {
                    "type": ["integer", "null"],
                    "minimum": 0,
                    "maximum": len(data["keys"]) - 1,
                }
                schema_values["items"]["items"]["items"][
                    field["category_slot"]
                ] = category
            except KeyError as err:
                show_error(err, 5)
    elif field["type"] == "array":
        schema_values["items"]["minItems"] = len(field["headers"])
        schema_values["items"]["maxItems"] = len(field["headers"])
        if "category_slot" in field:
            try:
                schema_values["items"]["items"] = [{} for x in field["headers"]]
                category = {
                    "type": ["integer", "null"],
                    "minimum": 0,
                    "maximum": len(data["keys"]) - 1,
                }
                schema_values["items"]["items"][field["category_slot"]] = category
            except KeyError as err:
                show_error(err, 5)
    return data_schema


def main(args):
    """Entry point for validator."""
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    schema_level = "basic" if args["--basic"] else "meta"
    schema_dir = os.path.join(os.path.dirname(script_dir), "data", "schema")
    schema = os.path.join(schema_dir, "%s.schema.json" % schema_level)
    file = args["DIRECTORY"]
    if os.path.isdir(file):
        file = os.path.join(file, "meta.json")
    if not os.path.isfile(file):
        if args["--example"]:
            print("Validating example file.")
            file = os.path.join(
                os.path.dirname(script_dir), "data", "example", "FXWY01", "meta.json"
            )
        else:
            show_error("%s not found" % file, 1)
    validate = fastjsonschema.compile(load_json_file(schema, True))

    meta = load_json_file(file)

    try:
        validate(meta)
    except fastjsonschema.exceptions.JsonSchemaException as err:
        show_error(err, 3)

    flattened_fields = flatten_fields(meta["fields"])
    try:
        meta["fields"] = flattened_fields
        validate(meta)
    except fastjsonschema.exceptions.JsonSchemaException as err:
        show_error(err, 4)

    validators = {}
    data_schemas = {}
    for type in ["array", "category", "identifier", "multiarray", "variable"]:
        type_schema = os.path.join(
            schema_dir, "subschemas", "%s.meta.schema.json" % type
        )
        validators[type] = fastjsonschema.compile(load_json_file(type_schema, True))
        data_schemas[type] = load_json_file(
            os.path.join(schema_dir, "subschemas", "%s.data.schema.json" % type), True,
        )

    for field in meta["fields"]:
        try:
            validators[field["type"]](field)
        except fastjsonschema.exceptions.JsonSchemaException as err:
            show_error(err, 4)
        data_file = file.replace("meta.json", "%s.json" % field["id"])
        data = load_json_file(data_file)
        data_schema = generate_data_schema(field, meta, data_schemas, data)
        data_validator = fastjsonschema.compile(data_schema)
        try:
            data_validator(data)
        except fastjsonschema.exceptions.JsonSchemaException as err:
            print(field["id"])
            show_error(err, 6)
    print("VALID")


def cli():
    """Entry point."""
    if len(sys.argv) == sys.argv.index(__name__.split(".")[-1]) + 1:
        args = docopt(__doc__, argv=[])
    else:
        args = docopt(__doc__, version=__version__)
    main(args)


if __name__ == "__main__":
    cli()
