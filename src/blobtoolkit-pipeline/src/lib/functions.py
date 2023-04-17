import pathlib


def reads_by_prefix(config):
    """Return read meta by prefix"""
    reads = {}
    if "reads" not in config:
        return {}

    for strategy in ("paired", "single"):
        if not strategy in config["reads"] or not config["reads"][strategy]:
            continue
        for entry in config["reads"][strategy]:
            if isinstance(entry, list):
                meta = {
                    "prefix": entry[0],
                    "platform": entry[1],
                    "strategy": strategy,
                    "base_count": entry[2],
                    "file": entry[3],
                }
                if len(entry) == 5:
                    meta["url"] = entry[4].split(";")
            else:
                meta = entry
            reads.update({meta["prefix"]: meta})
    return reads


def minimap_tuning(config, prefix):
    """Set minimap2 mapping parameter."""
    reads = reads_by_prefix(config)
    tunings = {
        "ILLUMINA": "sr",
        "OXFORD_NANOPORE": "map-ont",
        "PACBIO_SMRT": "map-pb",
        "OTHER": "sr",
    }
    return tunings.get(reads[prefix]["platform"], "sr")


def read_files(config, prefix):
    """Get read filenames."""
    reads = reads_by_prefix(config)
    return reads[prefix]["file"].split(";")


def seqtk_sample_input(config, prefix):
    """Generate seqtk command to subsamplereads if required."""
    meta = reads_by_prefix(config)[prefix]
    filenames = meta["file"].split(";")
    ratio = 1
    if "coverage" in config["reads"] and "max" in config["reads"]["coverage"]:
        base_count = meta.get("base_count", None)
        assembly_span = config["assembly"].get("span", None)
        if isinstance(base_count, int) and isinstance(assembly_span, int):
            try:
                ratio = assembly_span * config["reads"]["coverage"]["max"] / base_count
            except ZeroDivisionError:
                pass
    if ratio <= 0.95:
        command = " ".join(
            [
                "<(seqtk sample -s 100 %s %.2f)" % (filename, ratio)
                for filename in filenames
            ]
        )
    else:
        command = " ".join(filenames)
    return command


def diamond_db_name(config):
    """Generate filtered diamond database name."""
    name = "reference_proteomes"
    parts = ["diamond", name]
    return ".".join(parts)


def blobdir_name(config):
    """Generate blobdir name."""
    name = config["assembly"]["prefix"]
    if "revision" in config and config["revision"] > 0:
        name = "%s.%d" % (name, config["revision"])
    return name


def blobtools_cov_flag(config):
    """Generate --cov flag for blobtools add."""
    keys = reads_by_prefix(config).keys()
    if keys:
        return "--cov " + " --cov ".join(
            [
                "%s/%s.%s.bam=%s"
                % (minimap_path, config["assembly"]["prefix"], key, key)
                for key in keys
            ]
        )
    return ""


def gzipped_bed_cols(config):
    """Generate cut commands for gzipped bed files."""
    keys = reads_by_prefix(config).keys()
    if keys:
        return " ".join(
            [
                "<(echo %s_cov && gunzip -c %s.%s.regions.bed.gz | cut -f 4)"
                % (key, config["assembly"]["prefix"], key)
                for key in keys
            ]
        )
    return ""


def set_stats_chunk(config):
    """Set chunk size for calculating stats."""
    return config["settings"].get("stats_chunk", 1000)


def skip_windowmasker(config):
    """Set flag to skip windowmasker step."""
    return config["settings"].get("skip_windowmasker", False)


def set_stats_windows(config):
    """Set window sizes for averaging stats."""
    windows = config["settings"].get("stats_windows", 0.1)
    if not isinstance(windows, list):
        windows = [windows]
    if 1 not in windows:
        windows.append(1)
    return windows


def set_blast_chunk(config):
    """Set minimum chunk size for splitting long sequences."""
    return config["settings"].get("blast_chunk", 100000)


def set_blast_chunk_overlap(config):
    """Set overlap length for splitting long sequences."""
    return config["settings"].get("blast_overlap", 0)


def set_blast_max_chunks(config):
    """Set minimum chunk size for splitting long sequences."""
    return config["settings"].get("blast_max_chunks", 10)


def set_blast_min_length(config):
    """Set minimum sequence length for running blast searches."""
    return config["settings"].get("blast_min_length", 1000)


def taxid_flag(config, tool):
    """Set taxid flag for running blast searches."""
    taxid = config["taxon"].get("taxid", None)
    if taxid is None:
        return ""
    flags = {
        "blastn": "-negative_taxids",
        "diamond_blastp": "--taxon-exclude",
        "diamond_blastx": "--taxon-exclude",
    }
    return "%s %s" % (flags[tool], taxid)


def read_similarity_settings(config, group):
    """Read similarity settings for blast rules and outputs."""
    settings = {
        "evalue": 1.0e-10,
        "import_evalue": 1.0e-25,
        "max_target_seqs": 10,
        "name": "reference_proteomes",
    }
    if "defaults" in config["similarity"]:
        settings.update({**config["similarity"]["defaults"]})
    if group in config["similarity"]:
        settings.update({**config["similarity"][group]})
    return settings


def similarity_setting(config, group, value):
    """Get a single similarity setting value."""
    settings = read_similarity_settings(config, group)
    setting = settings.get(value, None)
    if setting is not None:
        return setting
    settings = {
        "evalue": 1.0e-25,
        "max_target_seqs": 10,
        "taxrule": "buscogenes",
        **settings,
    }
    if value.startswith("import_"):
        value = value.replace("import_", "")
    return settings[value]


def set_update_taxrule(config):
    taxrule = similarity_setting(config, "diamond_blastx", "taxrule")
    if taxrule == "buscoregions":
        return "--update-plot"
    return ""


def set_fields(config):
    params = ""
    if "fields" in config and isinstance(config["fields"], dict):
        for key, obj in config["fields"].items():
            if "file" in obj:
                if key == "synonyms":
                    prefix = obj.get("prefix", "")
                    if prefix:
                        prefix = "=%s" % prefix
                    params += " --synonyms %s%s" % (obj["file"], prefix)
                else:
                    params += " --text %s" % (obj["file"])
    if params:
        params += " --text-header"
    return params


def get_basal_lineages(config):
    """Get basal BUSCO lineages from config."""
    if "basal_lineages" in config["busco"]:
        return config["busco"]["basal_lineages"]
    basal = {"archaea_odb10", "bacteria_odb10", "eukaryota_odb10"}
    lineages = []
    if "lineages" in config["busco"]:
        for lineage in config["busco"]["lineages"]:
            if lineage in basal:
                lineages.push(lineage)
    return lineages


def set_view_timeout(config):
    """Set timeout duration for generating static images."""
    return config["settings"].get("view_timeout", 60)
