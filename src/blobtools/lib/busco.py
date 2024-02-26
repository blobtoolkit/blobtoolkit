#!/usr/bin/env python3
"""Parse BUSCO results into MultiArray Field."""

import re
from collections import defaultdict

from ..lib import file_io
from .field import MultiArray


def parse_busco(busco_file, identifiers):  # pylint: disable=too-many-locals
    """Parse BUSCO results into a MultiArray."""
    data = file_io.read_file(busco_file)
    if not data:
        print("WARNING file %s is empty" % busco_file)
        return None
    lines = data.split("\n")
    version = lines[0].split(":")[1].strip()
    desc = re.split(r":\s*|\(|\)\s*|,\s*", lines[1])
    meta = {
        "version": version,
        "set": desc[1].strip(),
        "count": int(desc[7].strip()),
        "file": busco_file,
    }
    version = int(version.split(".")[0])
    if version < 4:
        rows = [re.split("\t", line) for line in lines[5:]]
        meta["set"] = re.search(
            r"-l\s.*?\/*(\w+_odb\d+)\/", lines[2].split(":")[1].strip()
        )[1]
        columns = re.split(r"# |\t", lines[4])[1:]
        try:
            contig_index = columns.index("Contig")
        except ValueError:
            contig_index = columns.index("Sequence")
    else:
        rows = [re.split("\t", line) for line in lines[3:]]
        columns = re.split(r"# |\t", lines[2])[1:]
        contig_index = columns.index("Sequence")
    meta["field_id"] = "%s_busco" % meta["set"]
    busco_index = columns.index("Busco id")
    status_index = columns.index("Status")
    results = defaultdict(list)
    for row in rows:
        if len(row) > contig_index:
            if version < 4:
                contig = row[contig_index]
            else:
                contig = re.sub(r":\d+-\d+$", "", row[contig_index])
            results[contig].append([row[busco_index], row[status_index]])
    if not identifiers.validate_list(list(results.keys())):
        # try removing _\d+ suffix added by prokka-based busco
        res = {}
        for contig, values in results.items():
            ctg = re.sub(r"_\d+$", "", contig)
            if ctg not in res:
                res[ctg] = values
            else:
                res[ctg] += values
        results = res
        if not identifiers.validate_list(list(results.keys())):
            raise UserWarning(
                "Contig names in the Busco file did not match dataset identifiers."
            )
    values = [results[id] if id in results else [] for id in identifiers.values]
    busco_field = MultiArray(
        meta["field_id"],
        values=values,
        meta=meta,
        headers=("Busco id", "Status"),
        parents=["children"],
        category_slot=1,
    )
    return busco_field


def parse(files, **kwargs):
    """Parse all BUSCO files."""
    parsed = []
    for file in files:
        busco = parse_busco(file, identifiers=kwargs["dependencies"]["identifiers"])
        if busco is not None:
            parsed.append(busco)
    return parsed


def parent():
    """Set standard metadata for BUSCO."""
    busco = {"datatype": "mixed", "type": "array", "id": "busco", "name": "Busco"}
    return [busco]


def busco_score(values, total):
    """Calculate BUSCO score."""
    fragmented = set()
    buscos = defaultdict(int)
    for contig in values:
        for busco in contig:
            if busco[1] == "Fragmented":
                fragmented.add(busco[0])
            buscos[busco[0]] += 1
    present = len(buscos.keys())
    complete = present - len(fragmented)
    duplicated = len([value for value in buscos.values() if value > 1])
    scores = {
        "c": complete,
        "d": duplicated,
        "m": total - present,
        "f": len(fragmented),
        "t": total,
        "s": complete - duplicated,
    }
    string = "C:{:.1%}".format(scores["c"] / total)
    string += "[S:{:.1%},".format(scores["s"] / total)
    string += "D:{:.1%}],".format(scores["d"] / total)
    string += "F:{:.1%},".format(scores["f"] / total)
    string += "M:{:.1%},".format(scores["m"] / total)
    string += "n:{:d}".format(total)
    scores.update({"string": string})
    return scores


def summarise(indices, fields, **kwargs):  # pylint: disable=unused-argument
    """Summarise BUSCOs."""
    summary = {}
    for lineage in fields["lineages"]:
        values = fields[lineage].expand_values()
        values = [values[i] for i in indices]
        total = fields[lineage].meta["count"]
        lineage = lineage.replace("_busco", "")
        summary.update({lineage: busco_score(values, total)})
    return summary


def remove_from_meta(meta):
    """Delete all BUSCO fields."""
    field_ids = []
    if meta.has_field("busco"):
        field_ids = meta.remove_field("busco")
    return field_ids
