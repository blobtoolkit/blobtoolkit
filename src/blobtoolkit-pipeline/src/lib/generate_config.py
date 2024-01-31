#!/usr/bin/env python3

"""
Generate config files for BlobToolKit pipeline.

Usage:
  blobtoolkit-pipeline generate-config <ACCESSION> [--coverage 30] [--download]
    [--out /path/to/output/directory] [--db /path/to/database/directory]
    [--db-suffix STRING] [--reads STRING...] [--read-runs INT] [--api-key STRING]
    [--platforms STRING] [--datasets] [--protocol STRING]
    [--download-client STRING] [--retry-times INT]

Options:
  --coverage=INT            Maximum coverage for read mapping [default: 30]
  --download                Flag to download remote files [default: False]
  --out PATH                Path to output directory [default: .]
  --db PATH                 Path to database directory [default: .]
  --db-suffix STRING        Database version suffix (e.g. 2021_06)
  --reads STRING            Read accession to include.
  --read-runs INT           Maximum number of read runs [default: 3]
  --api-key STRING          NCBI api key for use with edirect
  --platforms STRING        priority order for sequencing platforms
                            [default: PACBIO_SMRT,ILLUMINA_XTEN,ILLUMINA,OXFORD_NANOPORE,OTHER]
  --protocol STRING         Download files using ftp(s) or http(s) [default: ftp]
  --download-client STRING  Client to use for downloads (aria2c|curl|datasets) [default: curl]
  --retry-times INT         Number of times to retry a failed download [default: 0]
"""

# pyright: reportMissingModuleSource=false

import gzip
import os
import re
import sys
from glob import glob
from operator import itemgetter
from pathlib import Path
from shutil import rmtree
from shutil import which
from subprocess import PIPE
from subprocess import Popen
from subprocess import run
from zipfile import ZipFile

import requests
import ujson
import yaml
from defusedxml import ElementTree as ET
from docopt import DocoptExit
from docopt import docopt
from tolkein import tofetch
from tolkein import tofile
from tolkein import tolog

LOGGER = tolog.logger(__name__)

GOAT_API = "https://goat.genomehubs.org/api/v2"
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
    lineages = [
        f'{BUSCO_SETS[obj["taxon_id"]]}_odb10'
        for obj in ancestors
        if obj["taxon_id"] in BUSCO_SETS
    ]
    lineages.extend(("bacteria_odb10", "archaea_odb10"))
    return lineages


def fetch_file_curl(url, filename, protocol="ftp", retry_times=0):
    """fetch a remote file using curl."""
    filepath = Path(filename)
    if filepath.is_file():
        LOGGER.info("File exists, not overwriting")
        return
    if protocol:
        url = url.replace("ftp:", f"{protocol}:").replace("ftps:", f"{protocol}:")
    iteration = 0
    complete = False
    while iteration <= retry_times:
        cmd = ["curl", "-L", "-s", "--connect-timeout", "30"]
        if iteration:
            LOGGER.info("Connection interrupted, resuming download")
            cmd.extend(["-C", "-"])
        cmd.extend(
            [
                "-o",
                f"{filepath.parent}/{filepath.name}",
                url,
            ]
        )
        proc = run(cmd)
        if proc.returncode == 0:
            complete = True
            break
        iteration += 1
    if not complete:
        LOGGER.error(f"Unable to fetch {filename}")
        sys.exit(1)


def fetch_file_aria(url, filename, protocol="ftp", retry_times=0, binary="aria2c"):
    """fetch a remote file using aria2."""
    filepath = Path(filename)
    if filepath.is_file():
        LOGGER.info("File exists, not overwriting")
        return
    if protocol:
        url = url.replace("ftp:", f"{protocol}:").replace("ftps:", f"{protocol}:")
    iteration = 0
    complete = False
    while iteration <= retry_times:
        cmd = [binary, "-q", "-d", filepath.parent]
        if iteration:
            LOGGER.info("Connection interrupted, resuming download")
            cmd.append("-c")
        cmd.extend(["-o", filepath.name, url])
        proc = run(cmd)
        if proc.returncode == 0:
            complete = True
            break
        iteration += 1
    if not complete:
        LOGGER.error(f"Unable to fetch {filename}")
        sys.exit(1)


def fetch_ncbi_datasets_assembly(filename, accession, retry_times):
    """Fetch a remote assembly using NCBI datasets."""
    filepath = Path(filename)
    if filepath.is_file():
        LOGGER.info("File exists, not overwriting")
        return
    dataset_zip = f"{filepath.parents[0]}/{accession}_dataset.zip"
    download_dehydrated = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--dehydrated",
        "--filename",
        dataset_zip,
    ]
    proc = run(download_dehydrated)
    if proc.returncode == 0:
        LOGGER.error(f"Unable to fetch {accession} with NCBI datasets")
        sys.exit(1)
    with ZipFile(dataset_zip, "r") as zh:
        zh.extractall(filepath.parents[0])
    rehydrate = ["datasets", "rehydrate", "--directory", filepath.parents[0]]
    iteration = 0
    complete = False
    while iteration <= retry_times:
        proc = run(rehydrate)
        if proc.returncode == 0:
            complete = True
            break
        iteration += 1
    if not complete:
        LOGGER.error(f"Unable to rehydrate {accession} with NCBI datasets")
        sys.exit(1)
    assembly_dir = f"{filepath.parents[0]}/ncbi_dataset/data/{accession}/"
    LOGGER.info(f"Writing compressed file to {filename}")
    with open(glob(f"{assembly_dir}/*_genomic.fna")[0], "rb") as ifh:
        with gzip.open(filename, "wb") as ofh:
            ofh.writelines(ifh)
    os.remove(dataset_zip)
    os.remove(f"{filepath.parents[0]}/README.md")
    rmtree(f"{filepath.parents[0]}/ncbi_dataset")


def fetch_assembly_url(accession, api_key=None):
    """
    Fetch an assembly url using edirect.
    """
    eutils_env = os.environ.copy()
    if api_key and api_key is not None:
        eutils_env["NCBI_API_KEY"] = api_key
    esearch = Popen(
        f"esearch -db assembly -query {accession}",
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
        if basenames := re.findall(GCA_NAME, url):
            for basename in basenames:
                if accession in basename:
                    return f"{url}/{basename}_genomic.fna.gz"
    return None


def fetch_assembly_fasta(
    url, filename, accession, protocol, retry_times, download_client
):
    """Save assembly fasta file to local disk."""
    if download_client == "datasets" and which(download_client) is not None:
        LOGGER.info(f"Fetching assembly FASTA to {filename} using NCBI datasets")
        fetch_ncbi_datasets_assembly(filename, accession, retry_times)

    elif download_client.startswith("aria") and which(download_client) is not None:
        LOGGER.info(f"Fetching assembly FASTA to {filename} using {download_client}")
        fetch_file_aria(url, filename, protocol, retry_times, download_client)
    else:
        LOGGER.info(f"Fetching assembly FASTA to {filename} using curl")
        fetch_file_curl(url, filename, protocol, retry_times)


def write_cat_file(data, filename):
    result = [key for key, obj in data.items()]
    arr = [result]
    for idx, value in enumerate(data[result[0]]["list"]):
        row = [value]
        row.extend(data[key]["list"][idx] for key in result[1:])
        arr.append(row)
    tofile.write_file(filename, arr)
    return result


def parse_assembly_report(filename, cat_filename, syn_filename):
    """Parse synonyms and assembly level into tsv files."""
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
    write_cat_file(cats, cat_filename)
    write_cat_file(names, syn_filename)


def fetch_assembly_report(
    url, filename, cat_filename, syn_filename, protocol, retry_times
):
    """Save assembly report file to local disk."""
    LOGGER.info(f"Fetching assembly report to {filename}")
    fetch_file_curl(url, filename, protocol, retry_times)
    parse_assembly_report(filename, cat_filename, syn_filename)


def fetch_assembly_meta_xml(accession):
    """
    Fetch assembly metadata xml from ENA.
    """
    return tofetch.fetch_url(f"{ENA_API}/xml/{accession}")


def deep_find_text(data, tags):
    """
    Find nested attributes in xml.

    Return attribute value.
    """
    for tag in tags:
        try:
            data = data.find(tag)
        except Exception:
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
    try:
        root = ET.fromstring(xml)
    except TypeError:
        LOGGER.error(f"Invalid accession: '{accession}'")
        sys.exit(1)
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
        meta["assembly"]["prefix"] = f"{wgs_prefix}{wgs_version.zfill(2)}"
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
        busco_lineage = f"{buscodir}/lineages/{lineage}"
        if not os.path.isdir(busco_lineage):
            lineages_to_fetch.append(lineage)
    if not lineages_to_fetch:
        return
    lineage_urls = {}
    LOGGER.info("Fetching BUSCO lineage directory listing")
    listing = tofetch.fetch_url(f"{BUSCO_URL}/")
    for entry in listing.split("\n"):
        parts = re.split(r"[\"\s]+", entry)
        if len(parts) == 8:
            busco_set = re.sub(r"\..+$", "", parts[2])
            lineage_urls[busco_set] = f"{BUSCO_URL}/{parts[2]}"
    for lineage in lineages_to_fetch:
        LOGGER.info(f"Fetching BUSCO lineage {lineage}")
        tofetch.fetch_tar(lineage_urls[lineage], buscodir)


def fetch_goat_data(taxon_id):
    """Fetch taxon metadata from GoaT."""
    LOGGER.info("Fetching taxon metadata")
    url = f"{GOAT_API}/record?recordId=taxon_id-{taxon_id}&result=taxon"
    result = tofetch.fetch_url(url)
    if result is None:
        LOGGER.error("Unable to fetch taxon metadata for '%s' from GoaT", taxon_id)
        sys.exit(1)
    data = ujson.loads(result)
    return data["records"][0]["record"]


def fetch_read_info(accession, per_platform):
    """Fetch read info for an accession."""
    portal = "https://www.ebi.ac.uk/ena/portal/api"
    url = f"{portal}/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_bytes,base_count,library_strategy,library_selection,library_layout,instrument_platform,experiment_title,fastq_ftp"

    si_suffix = {
        "k": 3,
        "M": 6,
        "G": 9,
        "T": 12,
        "P": 15,
        "E": 18,
    }

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
        for i in range(len(header)):
            value = fields[i]
            if header[i] == "instrument_platform" and platform != "ILLUMINA_XTEN":
                platform = fields[i]
            if (
                header[i] == "experiment_title"
                and value == "HiSeq X Ten paired end sequencing"
            ):
                platform = "ILLUMINA_XTEN"
            values[header[i]] = value
        try:
            values["base_count"] = int(values["base_count"])
        except ValueError:
            if "base_count" in values and values["base_count"].endswith("b"):
                values["base_count"] = values["base_count"].strip("b")
                suffix = values["base_count"][-1]
                try:
                    exponent = si_suffix[suffix]
                except KeyError:
                    exponent = 0
                values["base_count"] = int(values["base_count"].strip(suffix)) * pow(
                    10, exponent
                )
            else:
                values["base_count"] = 0
        except KeyError:
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
        runs = []
        if arr := per_platform[key]:
            runs.extend(iter(arr))
        runs = sorted(runs, key=itemgetter("base_count"), reverse=True)[
            : read_runs - len(sra)
        ]
        sra += runs
        if len(sra) >= read_runs:
            break
    return sra or None


def base_count(x):
    """Return number of bases or zero."""
    return int(x["base_count"][0] or 0) if isinstance(x["base_count"], list) else 0


def fetch_read_files(meta, protocol, retry_times):
    """Fetch sra reads."""
    files = meta["file"].split(";")
    for index, url in enumerate(meta["fastq_ftp"].split(";")):
        url = f"ftp://{url}"
        read_file = files[index]
        LOGGER.info("Fetching read file %s", read_file)
        fetch_file_curl(url, read_file, protocol, retry_times)


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
    for library in sra:
        strategy = library["library_layout"].lower()
        fastq_ftp = library["fastq_ftp"].split(";")[-2:]
        info = {
            "prefix": library["run_accession"],
            "platform": library["instrument_platform"],
            "base_count": library["base_count"],
            "file": ";".join(
                [re.sub(r"^.+/", f"{readdir}/", url) for url in fastq_ftp]
            ),
            "url": ";".join(fastq_ftp),
        }

        library["file"] = info["file"]
        meta["reads"][strategy].append(info)


def set_btk_version(meta):
    """Get curent version of hosted datasets from BTK API."""
    LOGGER.info("Checking current version on BTK public viewer")
    string = meta["assembly"]["prefix"]
    btk = f"https://blobtoolkit.genomehubs.org/api/v1/search/{string}"
    response = requests.get(btk)
    current = 0
    if response.ok:
        data = yaml.full_load(response.text)
        for asm in data:
            current = asm["version"] if "version" in asm else 1
    meta["revision"] = current
    meta["version"] = current + 1


def process_reads(opts, outdir, protocol, retry_times, meta, sra):
    if opts["--coverage"]:
        meta["reads"].update({"coverage": {"max": int(opts["--coverage"])}})
    readdir = f"{outdir}/reads"
    add_reads_to_meta(meta, sra, readdir)
    if opts["--download"]:
        os.makedirs(readdir, exist_ok=True)
        for library in sra:
            fetch_read_files(library, protocol, retry_times)


def process_busco(opts, buscodir, meta, taxon_meta):
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


def download_assembly(
    accession,
    outdir,
    protocol,
    download_client,
    retry_times,
    buscodir,
    assembly_url,
    assembly_file,
    assembly_report,
    syn_filename,
    cat_filename,
):
    os.makedirs(buscodir, exist_ok=True)
    os.makedirs(f"{outdir}/assembly", exist_ok=True)
    fetch_assembly_fasta(
        assembly_url,
        assembly_file,
        accession,
        protocol,
        retry_times,
        download_client,
    )
    report_url = assembly_url.replace("_genomic.fna.gz", "_assembly_report.txt")
    fetch_assembly_report(
        report_url,
        assembly_report,
        cat_filename,
        syn_filename,
        protocol,
        retry_times,
    )


def set_defaults(opts):
    accession = opts["<ACCESSION>"]
    outdir = opts["--out"]
    dbdir = opts["--db"]
    protocol = opts["--protocol"]
    download_client = opts["--download-client"]
    retry_times = int(opts["--retry-times"])
    buscodir = f"{dbdir}/busco"
    uniprotdir = f"{dbdir}/uniprot"
    ntdir = f"{dbdir}/nt"
    taxdumpdir = f"{dbdir}/taxdump"
    if opts["--db-suffix"]:
        buscodir += f'_{opts["--db-suffix"]}'
        ntdir += f'_{opts["--db-suffix"]}'
        uniprotdir += f'_{opts["--db-suffix"]}'
        taxdumpdir += f'_{opts["--db-suffix"]}'
    if not outdir.endswith(accession):
        outdir += f"/{accession}"
    return (
        accession,
        outdir,
        protocol,
        download_client,
        retry_times,
        buscodir,
        uniprotdir,
        ntdir,
        taxdumpdir,
    )


def main(rename=None):
    """Entry point."""
    docs = __doc__
    if rename is not None:
        docs = docs.replace("blobtoolkit-pipeline", rename)
    try:
        opts = docopt(docs)
    except DocoptExit as e:
        raise DocoptExit from e
    (
        accession,
        outdir,
        protocol,
        download_client,
        retry_times,
        buscodir,
        uniprotdir,
        ntdir,
        taxdumpdir,
    ) = set_defaults(opts)
    os.makedirs(outdir, exist_ok=True)
    meta = parse_assembly_meta(accession)
    assembly_url = fetch_assembly_url(accession, opts["--api-key"])
    if assembly_url is None:
        LOGGER.error("Unable to find assembly URL")
        sys.exit(1)
    assembly_file = f"{outdir}/assembly/{accession}.fasta.gz"
    meta["assembly"].update({"file": assembly_file, "url": assembly_url})
    assembly_report = f"{outdir}/assembly/{accession}.report.txt"
    syn_filename = f"{outdir}/assembly/{accession}.synonyms.tsv"
    cat_filename = f"{outdir}/assembly/{accession}.categories.tsv"
    meta["fields"] = {
        "synonyms": {"file": syn_filename, "prefix": "insdc"},
        "categories": {"file": cat_filename},
    }
    if opts["--download"]:
        download_assembly(
            accession,
            outdir,
            protocol,
            download_client,
            retry_times,
            buscodir,
            assembly_url,
            assembly_file,
            assembly_report,
            syn_filename,
            cat_filename,
        )
    taxon_meta = fetch_goat_data(meta["taxon"]["taxid"])
    add_taxon_to_meta(meta, taxon_meta)
    set_btk_version(meta)
    process_busco(opts, buscodir, meta, taxon_meta)
    read_accessions = []
    if meta["assembly"]["biosample"]:
        read_accessions = [meta["assembly"]["biosample"]]
    if opts["--reads"]:
        read_accessions += opts["--reads"]
    if sra := assembly_reads(
        read_accessions, int(opts["--read-runs"]), opts["--platforms"]
    ):
        process_reads(opts, outdir, protocol, retry_times, meta, sra)
    meta["similarity"]["blastn"]["path"] = ntdir
    meta["similarity"]["diamond_blastx"]["path"] = uniprotdir
    meta["similarity"]["diamond_blastp"]["path"] = uniprotdir
    meta["settings"]["taxdump"] = taxdumpdir
    tofile.write_file(f"{outdir}/config.yaml", meta)


if __name__ == "__main__":
    main()
