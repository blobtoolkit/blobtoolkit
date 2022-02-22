#!/usr/bin/env python3

"""
Generate config files for BlobToolKit pipeline.

Usage:
  btk pipeline generate-config <ACCESSION> [--coverage 30] [--download]
    [--out /path/to/output/directory] [--db /path/to/database/directory]
    [--db-suffix STRING] [--reads STRING...] [--read-runs INT] [--api-key STRING]
    [--platforms STRING]

Options:
  --coverage=INT         Maximum coverage for read mapping [default: 30]
  --download             Flag to download remote files [default: False]
  --out PATH             Path to output directory [default: .]
  --db PATH              Path to database directory [default: .]
  --db-suffix STRING     Database version suffix (e.g. 2021_06)
  --reads STRING         Read accession to include.
  --read-runs INT        Maximum number of read runs [default: 3]
  --api-key STRING       NCBI api key for use with edirect
  --platforms STRING     priority order for sequencing platforms
                         [default: PACBIO_SMRT,ILLUMINA_XTEN,ILLUMINA,OXFORD_NANOPORE,OTHER]
"""

import os
import re
import sys
from operator import itemgetter
from pathlib import Path
from subprocess import PIPE
from subprocess import Popen
from subprocess import run

import requests
import ujson
import yaml
from defusedxml import ElementTree as ET
from docopt import docopt
from tolkein import tofetch
from tolkein import tofile
from tolkein import tolog

LOGGER = tolog.logger(__name__)

GOAT_API = "https://goat.genomehubs.org/api/v0.0.1"
ENA_API = "https://www.ebi.ac.uk/ena/browser/api"
GCA_NAME = re.compile(r"GCA_.+")
BUSCO_URL = "https://busco-data.ezlab.org/v5/data/lineages"


def find_busco_lineages(ancestors):
    """Work out which BUSCO sets to run for a given taxid."""
    LOGGER.info("Identifying relevant BUSCO lineages")
    BUSCO_SETS = {
        "422676": "aconoidasida",
        "7898": "actinopterygii",
        "5338": "agaricales",
        "155619": "agaricomycetes",
        "33630": "alveolata",
        "5794": "apicomplexa",
        "6854": "arachnida",
        "6656": "arthropoda",
        "4890": "ascomycota",
        "8782": "aves",
        "5204": "basidiomycota",
        "68889": "boletales",
        "3699": "brassicales",
        "134362": "capnodiales",
        "33554": "carnivora",
        "91561": "cetartiodactyla",
        "34395": "chaetothyriales",
        "3041": "chlorophyta",
        "5796": "coccidia",
        "28738": "cyprinodontiformes",
        "7147": "diptera",
        "147541": "dothideomycetes",
        "3193": "embryophyta",
        "33392": "endopterygota",
        "314146": "euarchontoglires",
        "33682": "euglenozoa",
        "2759": "eukaryota",
        "5042": "eurotiales",
        "147545": "eurotiomycetes",
        "9347": "eutheria",
        "72025": "fabales",
        "4751": "fungi",
        "314147": "glires",
        "1028384": "glomerellales",
        "5178": "helotiales",
        "7524": "hemiptera",
        "7399": "hymenoptera",
        "5125": "hypocreales",
        "50557": "insecta",
        "314145": "laurasiatheria",
        "147548": "leotiomycetes",
        "7088": "lepidoptera",
        "4447": "liliopsida",
        "40674": "mammalia",
        "33208": "metazoa",
        "6029": "microsporidia",
        "6447": "mollusca",
        "4827": "mucorales",
        "1913637": "mucoromycota",
        "6231": "nematoda",
        "33183": "onygenales",
        "9126": "passeriformes",
        "5820": "plasmodium",
        "92860": "pleosporales",
        "38820": "poales",
        "5303": "polyporales",
        "9443": "primates",
        "4891": "saccharomycetes",
        "8457": "sauropsida",
        "4069": "solanales",
        "147550": "sordariomycetes",
        "33634": "stramenopiles",
        "32523": "tetrapoda",
        "155616": "tremellomycetes",
        "7742": "vertebrata",
        "33090": "viridiplantae",
        "71240": "eudicots",
    }
    # line = grep_line(lineage_file, taxid)
    lineages = []
    for obj in ancestors:
        if obj["taxon_id"] in BUSCO_SETS:
            lineages.append("%s_odb10" % BUSCO_SETS[obj["taxon_id"]])
    lineages.append("bacteria_odb10")
    lineages.append("archaea_odb10")
    return lineages


def fetch_file(url, filename):
    """fetch a remote file using aria2."""
    filepath = Path(filename)
    if filepath.is_file():
        LOGGER.info("File exists, not overwriting")
        return
    cmd = [
        "aria2c",
        "-q",
        "-d",
        filepath.parent,
        "-o",
        filepath.name,
        url,
    ]
    run(cmd)


def fetch_assembly_url(accession, api_key=None):
    """
    Fetch an assembly url using edirect.
    """
    eutils_env = os.environ.copy()
    if api_key and api_key is not None:
        eutils_env["NCBI_API_KEY"] = api_key
    esearch = Popen(
        "esearch -db assembly -query %s" % accession,
        stdout=PIPE,
        shell=True,
        env=eutils_env,
    )
    esummary = Popen(
        "esummary", stdin=esearch.stdout, stdout=PIPE, shell=True, env=eutils_env
    )
    xtract = Popen(
        "xtract -pattern DocumentSummary -element FtpPath_GenBank",
        stdin=esummary.stdout,
        stdout=PIPE,
        shell=True,
    )
    xtract_stdout = xtract.communicate()[0].decode("utf-8").strip()
    for url in xtract_stdout.split("\n"):
        basenames = re.findall(GCA_NAME, url)
        if basenames:
            for basename in basenames:
                if accession in basename:
                    return "%s/%s_genomic.fna.gz" % (url, basename)
    return None


def fetch_assembly_fasta(url, filename):
    """Save assembly fasta file to local disk."""
    LOGGER.info("Fetching assembly FASTA to %s" % filename)
    fetch_file(url, filename)


def parse_assembly_report(filename, cat_filename, syn_filename):
    """Parse synonyms and assembly level into tsv files."""
    synonyms = []
    categories = []
    cats = {
        "identifier": {"index": 4, "list": []},
        "assembly_role": {"index": 1, "list": []},
        "assembly_level": {"index": 3, "list": []},
        "assembly_unit": {"index": 7, "list": []},
    }
    names = {
        "identifier": {"index": 4, "list": []},
        "name": {"index": 0, "list": []},
        "assigned_name": {"index": 2, "list": []},
        "refseq_accession": {"index": 6, "list": []},
    }
    with tofile.open_file_handle(filename) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            row = line.rstrip().split("\t")
            for group in (cats, names):
                for obj in group.values():
                    value = row[obj["index"]]
                    obj["list"].append(value)
    header = []
    for key, obj in cats.items():
        if len(set(obj["list"])) > 1:
            header.append(key)
    categories.append(header)
    for idx, value in enumerate(cats[header[0]]["list"]):
        row = [value]
        for key in header[1:]:
            row.append(cats[key]["list"][idx])
        categories.append(row)
    tofile.write_file(cat_filename, categories)
    header = []
    for key, obj in names.items():
        if len(set(obj["list"])) > 1:
            header.append(key)
    synonyms.append(header)
    for idx, value in enumerate(names[header[0]]["list"]):
        row = [value]
        for key in header[1:]:
            row.append(names[key]["list"][idx])
        synonyms.append(row)
    tofile.write_file(syn_filename, synonyms)


def fetch_assembly_report(url, filename, cat_filename, syn_filename):
    """Save assembly report file to local disk."""
    LOGGER.info("Fetching assembly report to %s" % filename)
    fetch_file(url, filename)
    parse_assembly_report(filename, cat_filename, syn_filename)


def fetch_assembly_meta_xml(accession):
    """
    Fetch assembly metadata xml from ENA.
    """
    url = "%s/xml/%s" % (ENA_API, accession)
    xml = tofetch.fetch_url(url)
    return xml


def deep_find_text(data, tags):
    """
    Find nested attributes in xml.

    Return attribute value.
    """
    for tag in tags:
        try:
            data = data.find(tag)
        except:
            return ""
    return data.text


def parse_assembly_meta(accession):
    """Return dict of metadata values for an assembly."""
    LOGGER.info("Fetching assembly metadata")
    meta = {
        "assembly": {"accession": accession},
        "busco": {"lineages": []},
        "reads": {"paired": [], "single": []},
        "revision": 0,
        "settings": {
            "tmp": "/tmp",
            "blast_chunk": 100000,
            "blast_max_chunks": 10,
            "blast_overlap": 0,
            "blast_min_length": 1000,
            "stats_chunk": 1000,
            "stats_windows": [0.1, 0.01, 100000, 1000000],
        },
        "similarity": {
            "defaults": {
                "evalue": 1.0e-10,
                "import_evalue": 1.0e-25,
                "max_target_seqs": 10,
                "taxrule": "buscogenes",
            },
            "diamond_blastx": {"name": "reference_proteomes"},
            "diamond_blastp": {
                "name": "reference_proteomes",
                "import_max_target_seqs": 100000,
                "taxrule": "blastp=buscogenes",
            },
            "blastn": {"name": "nt"},
        },
        "taxon": {},
        "version": 1,
    }
    xml = fetch_assembly_meta_xml(accession)
    root = ET.fromstring(xml)
    asm = root.find("ASSEMBLY")
    meta["assembly"]["bioproject"] = deep_find_text(
        asm, ("STUDY_REF", "IDENTIFIERS", "PRIMARY_ID")
    )
    meta["assembly"]["biosample"] = deep_find_text(
        asm, ("SAMPLE_REF", "IDENTIFIERS", "PRIMARY_ID")
    )
    meta["taxon"]["taxid"] = deep_find_text(asm, ("TAXON", "TAXON_ID"))
    meta["taxon"]["name"] = deep_find_text(asm, ("TAXON", "SCIENTIFIC_NAME"))
    meta["assembly"]["level"] = asm.find("ASSEMBLY_LEVEL").text
    meta["assembly"]["alias"] = asm.attrib["alias"]
    wgs_prefix = deep_find_text(asm, ("WGS_SET", "PREFIX"))
    wgs_version = deep_find_text(asm, ("WGS_SET", "VERSION"))
    if wgs_prefix and wgs_version:
        meta["assembly"]["prefix"] = "%s%s" % (wgs_prefix, wgs_version.zfill(2))
    elif " " not in meta["assembly"]["alias"]:
        meta["assembly"]["prefix"] = meta["assembly"]["alias"].replace(".", "_")
    else:
        meta["assembly"]["prefix"] = meta["assembly"]["accession"].replace(".", "_")
    attributes = asm.find("ASSEMBLY_ATTRIBUTES")
    for attribute in attributes.findall("ASSEMBLY_ATTRIBUTE"):
        if attribute.find("TAG").text == "total-length":
            meta["assembly"]["span"] = int(attribute.find("VALUE").text)
        elif attribute.find("TAG").text == "scaffold-count":
            meta["assembly"]["scaffold-count"] = int(attribute.find("VALUE").text)
    return meta


def fetch_busco_lineages(busco_sets, buscodir):
    """Fetch busco lineages."""
    if not busco_sets:
        return
    lineages_to_fetch = []
    for lineage in busco_sets:
        busco_lineage = "%s/lineages/%s" % (buscodir, lineage)
        if not os.path.isdir(busco_lineage):
            lineages_to_fetch.append(lineage)
    if not lineages_to_fetch:
        return
    lineage_urls = {}
    LOGGER.info("Fetching BUSCO lineage directory listing")
    listing = tofetch.fetch_url("%s/" % BUSCO_URL)
    for entry in listing.split("\n"):
        parts = re.split(r"[\"\s]+", entry)
        if len(parts) == 8:
            busco_set = re.sub(r"\..+$", "", parts[2])
            lineage_urls.update({busco_set: "%s/%s" % (BUSCO_URL, parts[2])})
    for lineage in lineages_to_fetch:
        LOGGER.info("Fetching BUSCO lineage %s" % lineage)
        tofetch.fetch_tar(lineage_urls[lineage], buscodir)


def fetch_goat_data(taxon_id):
    """Fetch taxon metadata from GoaT."""
    LOGGER.info("Fetching taxon metadata")
    url = "%s/record?recordId=taxon_id-%s&result=taxon" % (GOAT_API, taxon_id)
    result = tofetch.fetch_url(url)
    if result is None:
        LOGGER.error("Unable to fetch taxon metadata for '%s' from GoaT", taxon_id)
        sys.exit(1)
    data = ujson.loads(result)
    return data["records"][0]["record"]


def fetch_read_info(accession, per_platform):
    """Fetch read info for an accession."""
    portal = "https://www.ebi.ac.uk/ena/portal/api"
    url = (
        "%s/filereport?accession=%s&result=read_run&fields=run_accession,fastq_bytes,base_count,library_strategy,library_selection,library_layout,instrument_platform,experiment_title,fastq_ftp"
        % (portal, accession)
    )
    data = tofetch.fetch_url(url)
    if data is None:
        return
    header = None
    for line in data.split("\n"):
        if not line or line == "":
            continue
        if header is None:
            header = line.split("\t")
            continue
        fields = line.split("\t")
        values = {}
        platform = "OTHER"
        for i in range(0, len(header)):
            value = fields[i]
            if header[i] == "instrument_platform" and platform != "ILLUMINA_XTEN":
                platform = fields[i]
            if header[i] == "experiment_title":
                if value == "HiSeq X Ten paired end sequencing":
                    platform = "ILLUMINA_XTEN"
            values.update({header[i]: value})
        if "base_count" in values:
            values["base_count"] = int(values["base_count"])
        else:
            values["base_count"] = 0
        try:
            per_platform[platform].append(values)
        except KeyError:
            per_platform["OTHER"].append(values)


def assembly_reads(accession, read_runs, platforms):
    """
    Query INSDC reads for an <accession> (or list of accessions).

    Return a dict of SRA accession, FASTQ ftp url, md5 and file size.
    """
    sra = []
    per_platform = {
        "PACBIO_SMRT": [],
        "ILLUMINA_XTEN": [],
        "ILLUMINA": [],
        "OXFORD_NANOPORE": [],
        "OTHER": [],
    }
    if not isinstance(accession, list):
        accession = [accession]
    for acc in accession:
        fetch_read_info(acc, per_platform)
    for key in platforms.split(","):
        arr = per_platform[key]
        runs = []
        if arr:
            for entry in arr:
                runs.append(entry)
        runs = sorted(runs, key=itemgetter("base_count"), reverse=True)[
            : read_runs - len(sra)
        ]
        sra += runs
        if len(sra) >= read_runs:
            break
    if sra:
        return sra
    return None


def base_count(x):
    """Return number of bases or zero."""
    if isinstance(x["base_count"], list):
        return int(x["base_count"][0] or 0)
    else:
        return 0


def fetch_read_files(meta):
    """Fetch sra reads."""
    files = meta["file"].split(";")
    for index, url in enumerate(meta["fastq_ftp"].split(";")):
        url = "ftp://%s" % url
        read_file = files[index]
        LOGGER.info("Fetching read file %s", read_file)
        fetch_file(url, read_file)


def add_taxon_to_meta(meta, taxon_meta):
    """Add taxon info to metadata."""
    LOGGER.info("Adding taxon metadata to assembly metadata")
    ranks = [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
        "superkingdom",
    ]
    for obj in taxon_meta["lineage"]:
        if obj["taxon_rank"] in ranks:
            meta["taxon"].update({obj["taxon_rank"]: obj["scientific_name"]})


def add_reads_to_meta(meta, sra, readdir):
    """Add read accessions to metadata."""
    LOGGER.info("Adding read accessions to assembly metadata")
    for index, library in enumerate(sra):
        strategy = library["library_layout"].lower()
        fastq_ftp = library["fastq_ftp"].split(";")[-2:]
        info = {
            "prefix": library["run_accession"],
            "platform": library["instrument_platform"],
            "base_count": library["base_count"],
            "file": ";".join(
                [re.sub(r"^.+/", "%s/" % readdir, url) for url in fastq_ftp]
            ),
            "url": ";".join(fastq_ftp),
        }
        library["file"] = info["file"]
        meta["reads"][strategy].append(info)


def set_btk_version(meta):
    """Get curent version of hosted datasets from BTK API."""
    LOGGER.info("Checking current version on BTK public viewer")
    string = meta["assembly"]["prefix"]
    btk = "https://blobtoolkit.genomehubs.org/api/v1/search/%s" % string
    response = requests.get(btk)
    current = 0
    if response.ok:
        data = yaml.full_load(response.text)
        for asm in data:
            if "version" in asm:
                current = asm["version"]
            else:
                current = 1
    meta["revision"] = current
    meta["version"] = current + 1


# def create_outdir(span, version=1, _lineage="all"):
#     """Create output directory."""
#     span = tolkein.tobin.readable_bin(span)
#     name = "%s/v%s/%s" % (OUTDIR, str(version), span)
#     os.makedirs("%s/" % name, exist_ok=True)
#     return name


def main():
    """Entry point."""
    opts = docopt(__doc__)
    accession = opts["<ACCESSION>"]
    outdir = opts["--out"]
    dbdir = opts["--db"]
    buscodir = "%s/busco" % dbdir
    uniprotdir = "%s/uniprot" % dbdir
    ntdir = "%s/nt" % dbdir
    taxdumpdir = "%s/taxdump" % dbdir
    if opts["--db-suffix"]:
        buscodir += "_%s" % opts["--db-suffix"]
        ntdir += "_%s" % opts["--db-suffix"]
        uniprotdir += "_%s" % opts["--db-suffix"]
        taxdumpdir += "_%s" % opts["--db-suffix"]
    if not outdir.endswith(accession):
        outdir += "/%s" % accession
    os.makedirs(outdir, exist_ok=True)
    meta = parse_assembly_meta(accession)
    assembly_url = fetch_assembly_url(accession, opts["--api-key"])
    if assembly_url is None:
        LOGGER.error("Unable to find assembly URL")
        sys.exit(1)
    assembly_file = "%s/assembly/%s.fasta.gz" % (outdir, accession)
    meta["assembly"].update({"file": assembly_file, "url": assembly_url})
    assembly_report = "%s/assembly/%s.report.txt" % (outdir, accession)
    syn_filename = "%s/assembly/%s.synonyms.tsv" % (outdir, accession)
    cat_filename = "%s/assembly/%s.categories.tsv" % (outdir, accession)
    meta["fields"] = {
        "synonyms": {"file": syn_filename, "prefix": "insdc"},
        "categories": {"file": cat_filename},
    }
    if opts["--download"]:
        os.makedirs(buscodir, exist_ok=True)
        os.makedirs("%s/assembly" % outdir, exist_ok=True)
        fetch_assembly_fasta(assembly_url, assembly_file)
        report_url = assembly_url.replace("_genomic.fna.gz", "_assembly_report.txt")
        fetch_assembly_report(report_url, assembly_report, cat_filename, syn_filename)
    taxon_meta = fetch_goat_data(meta["taxon"]["taxid"])
    add_taxon_to_meta(meta, taxon_meta)
    set_btk_version(meta)
    busco_sets = find_busco_lineages(taxon_meta["lineage"])
    if busco_sets:
        meta["busco"].update(
            {
                "download_dir": buscodir,
                "lineages": busco_sets,
                "basal_lineages": [
                    "eukaryota_odb10",
                    "bacteria_odb10",
                    "archaea_odb10",
                ],
            }
        )
    if opts["--download"]:
        fetch_busco_lineages(busco_sets, buscodir)
    read_accessions = []
    if meta["assembly"]["biosample"]:
        read_accessions = [meta["assembly"]["biosample"]]
    if opts["--reads"]:
        read_accessions += opts["--reads"]
    sra = assembly_reads(read_accessions, int(opts["--read-runs"]), opts["--platforms"])
    if sra:
        if opts["--coverage"]:
            meta["reads"].update({"coverage": {"max": int(opts["--coverage"])}})
        readdir = "%s/reads" % outdir
        add_reads_to_meta(meta, sra, readdir)
        if opts["--download"]:
            os.makedirs(readdir, exist_ok=True)
            for library in sra:
                fetch_read_files(library)
    meta["similarity"]["blastn"].update({"path": ntdir})
    meta["similarity"]["diamond_blastx"].update({"path": uniprotdir})
    meta["similarity"]["diamond_blastp"].update({"path": uniprotdir})
    meta["settings"]["taxdump"] = taxdumpdir
    tofile.write_file("%s/config.yaml" % outdir, meta)


if __name__ == "__main__":
    main()
