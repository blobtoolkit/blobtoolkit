#!/usr/bin/env python3

# pylint: disable=too-many-locals,too-many-branches,unused-argument,unused-variable
# pylint: disable=too-many-arguments,too-many-statements

"""Parse BLAST results into MultiArray Field."""

import math
import re
from collections import Counter
from collections import defaultdict
from itertools import groupby
from operator import itemgetter

from ..lib import file_io
from .field import Category
from .field import MultiArray
from .field import Variable


def parse_blast(blast_file, cols, results=None, index=0, evalue=1, bitscore=1):
    """Parse file into dict of lists."""
    if results is None:
        results = defaultdict(list)
    bitscores = {}
    blastp = {}

    for line in file_io.stream_file(blast_file):
        row = line.rstrip().split("\t")
        score = float(row[cols["bitscore"]])
        if score < bitscore:
            continue
        if len(row) == 4:
            cols["sseqid"] = 3
        else:
            if evalue < float(row[cols["evalue"]]):
                continue
        # allow for mis-specified columns following documentation bug
        if "sstart" in cols and "qstart" not in cols:
            cols["qstart"] = cols["sstart"]
        if "send" in cols and "qend" not in cols:
            cols["qend"] = cols["send"]
        seq_id, *offset = row[cols["qseqid"]].split("_-_")
        offset = int(offset[0]) if offset else 0
        query = row[cols["qseqid"]]
        if ":" in query and "=" in query:
            # parse blastp
            parts = query.split("=")
            if query in bitscores and score <= bitscores[query]:
                continue
            if len(parts) == 3 and parts[2] == "fragmented":
                continue
            bitscores[query] = score
            seq_id, start, end = re.split(r"[:-]", parts[0])
            hit = {
                "subject": row[cols["sseqid"]],
                "score": score,
                "start": int(start),
                "end": int(end),
                "file": index,
                "title": parts[1],
            }
        else:
            # parse blastx/blastn
            try:
                hit = {
                    "subject": row[cols["sseqid"]],
                    "score": score,
                    "start": int(row[cols["qstart"]]) + offset,
                    "end": int(row[cols["qend"]]) + offset,
                    "file": index,
                }
            except IndexError:
                # parse file without positions
                hit = {
                    "subject": row[cols["sseqid"]],
                    "score": score,
                    "start": None,
                    "end": None,
                    "file": index,
                }
        try:
            taxid = row[cols["staxids"]]
            try:
                taxid, *rest = taxid.split(";")
            except ValueError:
                # no taxid for this row
                pass
            hit.update({"taxid": int(taxid)})
        except ValueError:
            # no taxid in file
            hit.update({"taxid": 0})
        if bitscores:
            blastp[query] = hit
        else:
            results[seq_id].append(hit)
    if bitscores:
        for query, hit in blastp.items():
            seq_id, rest = query.split(":")
            results[seq_id].append(hit)
    return results


def chunk_size(value):
    """Calculate nice value for chunk size."""
    mag = math.floor(math.log10(value))
    first = int(str(value)[:2]) + 1
    size = first * pow(10, mag - 1)
    return size


def set_windows(meta, taxrule):
    """Set window sizes from dataset metadata."""
    windows = []
    chunk = 0.1
    try:
        if "bestsum" in taxrule:
            chunk = 1
        elif "blast_max_chunks" in meta.settings:
            chunk = 1 / int(meta.settings["blast_max_chunks"])
        if "stats_windows" in meta.settings:
            for x in meta.settings["stats_windows"]:
                window = {
                    "key": ("%f" % x).rstrip("0").rstrip("."),
                    "window": True,
                    "value": float(x),
                    "chunk": False,
                }
                if window["key"] == "0.1":
                    window["title"] = "windows"
                else:
                    window["title"] = "windows_%s" % window["key"]
                if chunk is not None and float(x) == chunk:
                    window.update({"chunk": True})
                    chunk = None
                windows.append(window)
    except AttributeError:
        pass
    if chunk is not None:
        windows.append(
            {
                "key": ("%f" % chunk).rstrip("0").rstrip("."),
                "window": False,
                "value": chunk,
                "chunk": True,
            }
        )
    return windows


def bin_hits(
    blast,
    meta,
    identifiers,
    lengths,
    taxrule,
    hit_count,
):
    """Place hits into bins for each window size."""
    windows = set_windows(meta, taxrule)
    bins = defaultdict(dict)
    get_item = itemgetter("start")
    get_score = itemgetter("score")
    for idx, identifier in enumerate(identifiers.values):
        length = lengths[idx]
        hits = sorted(blast[identifier], key=get_item)
        for window in windows:
            if window["value"] > 1:
                chunk = window["value"]
            elif window["value"] == 1:
                chunk = length
            else:
                chunk = round(length * window["value"] / 1000 + 0.5) * 1000
            if length < 1000000:
                if not window["chunk"]:
                    continue
                chunk = 100000
            bin = {**window}
            groups = {
                k: list(g)
                for k, g in groupby(hits, lambda x: get_item(x) // chunk * chunk)
            }
            current = 0
            bin["bins"] = []
            while current < length:
                if current in groups:
                    bin["bins"].append(
                        sorted(groups[current], key=get_score, reverse=True)[:hit_count]
                    )
                else:
                    bin["bins"].append([])
                current += chunk
            bins[identifier][window["key"]] = bin
    return bins


def initialise_values(bins, prefix, results, ranks):
    """Initialise results and values ready to apply taxrule."""
    windows = [
        obj["title"]
        for obj in bins[list(bins.keys())[0]].values()
        if "title" in obj and obj["window"]
    ]
    if results is None:
        blank = [None] * len(bins.keys())
        results = []
        for rank in ranks:
            result = {
                "field_id": "%s_%s" % (prefix, rank),
                "values": blank[:],
                "data": {
                    "cindex": blank[:],
                    "score": blank[:],
                    "positions": blank[:],
                    "hits": blank[:],
                },
            }
            for window in windows:
                result["data"].update({window: blank[:]})
            results.append(result)
    values = []
    for rank in ranks:
        value = {
            "category": defaultdict(str),
            "cindex": defaultdict(int),
            "score": defaultdict(float),
            "positions": defaultdict(list),
            "categories": defaultdict(list),
            "scores": defaultdict(list),
            "hits": defaultdict(list),
        }
        for window in windows:
            value.update({window: defaultdict(list)})
        values.append(value)
    return results, values


def apply_taxrule_to_bin(seqid, values, bin, window, taxdump, ranks):
    """Apply taxrule to a single bin."""
    for index, rank in enumerate(ranks):
        cat_scores = defaultdict(float)
        for hit in bin:
            try:
                category = taxdump.ancestors[hit["taxid"]][rank]
            except KeyError:
                category = 0
            if category > 0:
                category = taxdump.names[category]
            elif category < 0:
                category = "%s-undef" % taxdump.names[-category]
            else:
                category = "undef"
            # category = "no-hit"
            if category != "undef":
                cat_scores[category] += hit["score"]
            if index == 0 and window["chunk"]:
                try:
                    values[index]["hits"][seqid].append(
                        [
                            hit["taxid"],
                            hit["start"],
                            hit["end"],
                            hit["score"],
                            hit["subject"],
                            hit["file"],
                            hit.get("title", None),
                        ]
                    )
                except KeyError:
                    pass
            if window["chunk"]:
                values[index]["positions"][seqid].append([category])
        if cat_scores:
            category = max(cat_scores, key=cat_scores.get)
            score = cat_scores[category]
        else:
            category = None
            score = 0
        if "title" in window and window["title"] in values[index]:
            values[index][window["title"]][seqid].append([category])
        if window["chunk"]:
            values[index]["categories"][seqid].append([category])
            values[index]["scores"][seqid].append(score)


def apply_taxrule_across_bins(seqid, values, window, ranks):
    """Find most common bin category."""
    if window["chunk"]:
        for index, rank in enumerate(ranks):
            scores = defaultdict(float)
            counts = defaultdict(int)
            for idx, category in enumerate(values[index]["categories"][seqid]):
                if category and category[0] is not None:
                    counts[category[0]] += 1
                    scores[category[0]] += values[index]["scores"][seqid][idx]
            max_score = 0
            max_count = 0
            top_cat = False
            for category, count in counts.items():
                if count > max_count:
                    max_count = count
                    max_score = scores[category]
                    top_cat = category
                elif count == max_count:
                    score = scores[category]
                    if score > max_score:
                        max_count = count
                        max_score = score
                        top_cat = category
            if top_cat is not False:
                values[index]["category"][seqid] = top_cat
                values[index]["cindex"][seqid] = len(counts.keys()) - 1
                values[index]["score"][seqid] = max_score


def add_values_to_results(results, values, ranks, identifiers):
    """Add new values to result set."""
    for index, rank in enumerate(ranks):
        if not identifiers.validate_list(list(values[index]["category"].keys())):
            raise UserWarning(
                "Contig names in the hits file do not match dataset identifiers."
            )
        for i, seq_id in enumerate(identifiers.values):
            if results[index]["data"]["score"][i] in (0, None):
                if seq_id in values[index]["category"]:
                    results[index]["values"][i] = values[index]["category"][seq_id]
                    results[index]["data"]["score"][i] = values[index]["score"][seq_id]
                    results[index]["data"]["cindex"][i] = values[index]["cindex"][
                        seq_id
                    ]
                    results[index]["data"]["positions"][i] = values[index]["positions"][
                        seq_id
                    ]
                    if index == 0:
                        results[index]["data"]["hits"][i] = values[index]["hits"][
                            seq_id
                        ]
                else:
                    results[index]["values"][i] = "no-hit"
                    results[index]["data"]["score"][i] = 0
                    results[index]["data"]["cindex"][i] = 0
                    results[index]["data"]["positions"][i] = []
                    if index == 0:
                        results[index]["data"]["hits"][i] = []
                for key in values[index].keys():
                    if key.startswith("windows"):
                        results[index]["data"][key][i] = values[index][key][seq_id]


def apply_taxrule(bins, taxdump, taxrule, prefix, results, identifiers):
    """Apply taxrule to binned BLAST results."""
    ranks = taxdump.list_ranks()
    # ranks = ["superkingdom"]
    results, values = initialise_values(bins, prefix, results, ranks)
    for seqid, windows in bins.items():
        for key, window in windows.items():
            for bin in window["bins"]:
                apply_taxrule_to_bin(seqid, values, bin, window, taxdump, ranks)
            apply_taxrule_across_bins(seqid, values, window, ranks)
    add_values_to_results(results, values, ranks, identifiers)
    return results


def create_fields(results, taxrule, files, fields=None):
    """Store BLAST results as Fields."""
    if fields is None:
        fields = []
    hits_id = "%s_%s" % (taxrule, "positions")
    fields.append(
        MultiArray(
            hits_id,
            values=results[0]["data"]["hits"],
            meta={
                "field_id": hits_id,
                "name": hits_id,
                "type": "multiarray",
                "datatype": "mixed",
                "preload": False,
                "active": False,
                "files": files,
            },
            parents=["children", {"id": taxrule}, "children"],
            category_slot=None,
            headers=["taxid", "start", "end", "score", "subject", "index", "title"],
        )
    )
    for result in results:
        main = Category(
            result["field_id"],
            values=result["values"],
            meta={"field_id": result["field_id"], "name": result["field_id"]},
            parents=["children", {"id": taxrule}, "children"],
        )
        fields.append(main)
        parents = [
            "children",
            {"id": taxrule},
            "children",
            {"id": result["field_id"]},
            "data",
        ]
        field_id = "%s_%s" % (result["field_id"], "cindex")
        fields.append(
            Variable(
                field_id,
                values=result["data"]["cindex"],
                meta={
                    "scale": "scaleLinear",
                    "field_id": field_id,
                    "name": field_id,
                    "datatype": "integer",
                    "range": [
                        min(result["data"]["cindex"]),
                        max(result["data"]["cindex"]),
                    ],
                    "preload": False,
                    "active": False,
                },
                parents=parents,
            )
        )
        field_id = "%s_%s" % (result["field_id"], "score")
        _min = min(result["data"]["score"])
        fields.append(
            Variable(
                field_id,
                values=result["data"]["score"],
                meta={
                    "scale": "scaleLog",
                    "field_id": field_id,
                    "name": field_id,
                    "clamp": 1 if _min == 0 else False,
                    "datatype": "float",
                    "range": [_min, max(result["data"]["score"])],
                    "preload": False,
                    "active": False,
                },
                parents=parents,
            )
        )
        subfield = "positions"
        field_id = "%s_%s" % (result["field_id"], subfield)
        if len(result["data"][subfield]) > 1:
            headers = ["name"]
        else:
            headers = ["name"]
        fields.append(
            MultiArray(
                field_id,
                values=result["data"][subfield],
                fixed_keys=main.keys,
                meta={
                    "field_id": field_id,
                    "name": field_id,
                    "type": "multiarray",
                    "datatype": "string",
                    "preload": False,
                    "active": False,
                    "linked_field": hits_id,
                },
                parents=parents,
                category_slot=0,
                headers=headers,
            )
        )
        for subfield in result["data"].keys():
            if subfield.startswith("windows"):
                field_id = "%s_%s" % (result["field_id"], subfield)
                if len(result["data"][subfield]) > 1:
                    headers = ["name"]
                else:
                    headers = ["name"]
                fields.append(
                    MultiArray(
                        field_id,
                        values=result["data"][subfield],
                        fixed_keys=main.keys,
                        meta={
                            "field_id": field_id,
                            "name": field_id,
                            "type": "array",
                            "datatype": "string",
                            "preload": False,
                            "active": False,
                        },
                        parents=parents,
                        category_slot=0,
                        headers=headers,
                    )
                )

    return fields


def parse(files, **kwargs):
    """Parse BLAST results into Fields."""
    blast = None
    fields = []
    identifiers = kwargs["dependencies"]["identifiers"]
    lengths = kwargs["dependencies"]["length"].values
    try:
        taxrule, prefix = kwargs["--taxrule"].split("=")
    except ValueError:
        taxrule = kwargs["--taxrule"]
        prefix = taxrule
    cols = {}
    columns = kwargs["--hits-cols"].split(",")
    for column in columns:
        try:
            index, name = column.split("=")
            cols[name] = int(index) - 1
        except ValueError:
            exit("ERROR: --hits-cols contains an invalid value.")
    if taxrule.endswith("order"):
        results = None
        for index, file in enumerate(files):
            blast = parse_blast(
                file,
                cols,
                None,
                index,
                float(kwargs["--evalue"]),
                float(kwargs["--bitscore"]),
            )
            bins = bin_hits(
                blast,
                kwargs["meta"],
                identifiers,
                lengths,
                taxrule,
                int(kwargs["--hit-count"]),
            )
            results = apply_taxrule(
                bins, kwargs["taxdump"], taxrule, prefix, results, identifiers
            )
        fields = create_fields(results, prefix, files)
    else:
        results = None
        previous = None
        for index, file in enumerate(files):
            blast = parse_blast(
                file,
                cols,
                previous,
                index,
                float(kwargs["--evalue"]),
                float(kwargs["--bitscore"]),
            )
            previous = blast
        bins = bin_hits(
            blast,
            kwargs["meta"],
            identifiers,
            lengths,
            taxrule,
            int(kwargs["--hit-count"]),
        )
        results = apply_taxrule(
            bins, kwargs["taxdump"], taxrule, prefix, results, identifiers
        )
        fields = create_fields(results, prefix, files)
    update_plot = kwargs.get("--update-plot", False)
    if update_plot or "cat" not in kwargs["meta"].plot:
        kwargs["meta"].plot.update({"cat": "%s_phylum" % prefix})
    return fields


def parent():
    """Set standard metadata for BLAST."""
    blast = {
        "datatype": "string",
        "type": "category",
        "id": "taxonomy",
        "name": "Taxonomy",
    }
    return [blast]


def summarise(indices, fields, **kwargs):
    """Summarise assembly sequence stats."""
    taxonomy = kwargs["stats"]["taxonomy"]
    summary = {}
    lengths = [fields["length"].values[i] for i in indices]
    gcs = [fields["gc"].values[i] for i in indices]
    covs = [fields["cov"].values[i] for i in indices] if "cov" in fields else []
    summary.update({"total": length_stats(lengths, gcs, covs)})
    hits = list(
        map(
            lambda x: fields["hits"].keys[x],
            [fields["hits"].values[i] for i in indices],
        )
    )
    counts = Counter(hits)
    index = 0
    limit = 9
    other = []
    other_gcs = []
    other_covs = []
    for key in sorted(counts.items(), key=lambda x: x[1], reverse=True):
        taxon = key[0]
        subset = [i for i, j in enumerate(hits) if j == taxon]
        lengths = [fields["length"].values[indices[i]] for i in subset]
        gcs = [fields["gc"].values[indices[i]] for i in subset]
        covs = (
            [fields["cov"].values[indices[i]] for i in subset]
            if "cov" in fields
            else []
        )
        if len(counts.keys()) > limit and index >= limit:
            other += lengths
            other_gcs += gcs
            other_covs += covs
            if "target" in taxonomy and taxon == taxonomy["target"]:
                summary.update({"target": length_stats(lengths, gcs, covs)})
        else:
            summary.update({taxon: length_stats(lengths, gcs, covs)})
        index += 1
    if other:
        summary.update({"other": length_stats(other, other_gcs, other_covs)})
    return summary


def weighted_mean_sd(values, weights, log=False):
    """Calculate weighted mean, median and standard deviation."""
    weighted = []
    total_weight = 0
    total = 0
    for weight, value in zip(weights, values):
        if log:
            value = math.log10(value) + 10
        weighted.append({"weight": weight, "value": value})
        total_weight += weight
        total += value * weight
    mid = total_weight / 2
    running_total = 0
    median = 0
    for entry in sorted(weighted, key=lambda i: i["value"]):
        running_total += entry["weight"]
        if running_total > mid:
            median = entry["value"]
            break
    mean = total / total_weight
    if median == 0:
        median = mean
    variance = (
        sum([entry["weight"] * ((entry["value"] - mean) ** 2) for entry in weighted])
        / total_weight
    )
    st_dev = variance ** 0.5
    upper = mean + 2 * st_dev
    lower = mean - 2 * st_dev
    if log:
        upper = 10 ** (upper - 10)
        lower = 10 ** (lower - 10)
        mean = 10 ** (mean - 10)
        median = 10 ** (median - 10)
    return mean, median, st_dev, upper, lower


def length_stats(all_lengths, all_gcs, all_covs):
    """Calculate stats for a set of sequence lengths."""
    span = sum(all_lengths)
    count = len(all_lengths)
    lengths = []
    gcs = []
    covs = []
    if all_covs:
        for length, gc_value, cov in zip(all_lengths, all_gcs, all_covs):
            if cov >= 0.01:
                lengths.append(length)
                gcs.append(gc_value)
                covs.append(cov)
    else:
        lengths = all_lengths
        gcs = all_gcs
    stats = {"span": span, "count": count}
    if gcs:
        gc_mean, gc_median, gc_dev, gc_upper, gc_lower = weighted_mean_sd(gcs, lengths)
        stats.update(
            {
                "gc": [
                    float("%.4f" % gc_mean),
                    float("%.4f" % gc_median),
                    float("%.4f" % gc_lower),
                    float("%.4f" % gc_upper),
                    float("%.4f" % min(gcs)),
                    float("%.4f" % max(gcs)),
                ]
            }
        )
    if covs:
        cov_mean, cov_median, cov_dev, cov_upper, cov_lower = weighted_mean_sd(
            covs, lengths, log=True
        )
        stats.update(
            {
                "cov": [
                    float("%.4f" % cov_mean),
                    float("%.4f" % cov_median),
                    float("%.4f" % cov_lower),
                    float("%.4f" % cov_upper),
                    float("%.4f" % min(covs)),
                    float("%.4f" % max(covs)),
                ]
            }
        )
    n50 = span * 0.5
    n90 = span * 0.9
    all_lengths.sort(reverse=True)
    nlength = 0
    ncount = 0
    for length in all_lengths:
        ncount += 1
        nlength += length
        if "n50" not in stats and nlength > n50:
            stats.update({"n50": length, "l50": ncount})
        if "n90" not in stats and nlength > n90:
            stats.update({"n90": length, "l90": ncount})
    if "n50" not in stats:
        stats.update({"n50": all_lengths[-1], "l50": ncount})
    if "n90" not in stats:
        stats.update({"n90": all_lengths[-1], "l90": ncount})
    return stats


def remove_from_meta(meta):
    """Delete all hits fields."""
    field_ids = []
    if meta.has_field("taxonomy"):
        field_ids += meta.remove_field("taxonomy")
        meta.plot.pop("cat", None)
    return field_ids
