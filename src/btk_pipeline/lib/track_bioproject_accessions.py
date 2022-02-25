#!/usr/bin/env python3

"""
Track new data releases under a bioproject.

Usage:
  track_bioproject_accessions.py <BIOPROJECT> [--out /path/to/output/directory]

Options:
  --out PATH             Path to output directory [default: .]
"""

import os
from time import sleep

from docopt import docopt
from tolkein import tofetch
from tolkein import tofile
from tolkein import tolog

LOGGER = tolog.logger(__name__)

ENA_API = "https://www.ebi.ac.uk/ena/portal/api"


def fetch_bioproject_children(
    bioproject, *, projects=None, children=None, new_projects=None
):
    """
    Fetch children of a bioproject.
    """
    if projects is None:
        projects = []
    if new_projects is None:
        new_projects = {}
    if children is None:
        children = []
    url = "%s/search?result=study&query=parent_study%%3D%%22%s%%22&format=tsv" % (
        ENA_API,
        bioproject,
    )
    result = tofetch.fetch_url(url)
    if result and result is not None:
        for line in result.split("\n")[1:]:
            if line and "\t" in line:
                child_accession, description = line.split("\t")
                if child_accession in projects:
                    continue
                if "genome assembly" in description:
                    if "alternate haplotype" not in description:
                        children.append(child_accession)
                        new_projects.update({child_accession: bioproject})
                else:
                    sleep(0.5)
                    LOGGER.info(
                        "Fetching nested accessions under bioproject %s"
                        % child_accession
                    )
                    fetch_bioproject_children(
                        child_accession,
                        projects=projects,
                        children=children,
                        new_projects=new_projects,
                    )
    return children, new_projects


def fetch_accession(bioproject):
    """Fetch a GCA accession for a bioproject."""
    LOGGER.info("Fetching GCA accession for bioproject %s" % bioproject)
    url = "%s/search?result=assembly&query=study_accession%%3D%%22%s%%22&fields=accession%%2Cversion&format=tsv" % (
        ENA_API,
        bioproject,
    )
    result = tofetch.fetch_url(url)
    accession = None
    if result and result is not None:
        line = result.split("\n")[1]
        if line and "\t" in line:
            accession, version = line.split("\t")
            accession += ".%s" % version
    return accession


def fetch_bioproject_accessions(bioproject, *, projects=None):
    """Save assembly fasta file to local disk."""
    if projects is None:
        projects = []
    LOGGER.info("Fetching accessions under bioproject %s" % bioproject)
    accessions = []
    children, new_projects = fetch_bioproject_children(bioproject, projects=projects)
    project_list = set()
    for child_bioproject in children:
        sleep(0.5)
        accession = fetch_accession(child_bioproject)
        if accession is not None:
            accessions.append(accession)
            # project_list.add(child_bioproject)
            project_list.add(new_projects[child_bioproject])
    return accessions, project_list


if __name__ == "__main__":
    opts = docopt(__doc__)
    bioproject = opts["<BIOPROJECT>"]
    outdir = opts["--out"]
    os.makedirs(outdir, exist_ok=True)
    projects_file = "%s/bioprojects.processed" % outdir
    targets_file = "%s/accessions.todo" % outdir
    active_file = "%s/accessions.active" % outdir
    processed_file = "%s/accessions.processed" % outdir
    projects = tofile.read_file(projects_file)

    accessions, new_projects = fetch_bioproject_accessions(
        bioproject, projects=projects
    )
    with open(targets_file, "a+") as tfh:
        tfh.writelines([accession + "\n" for accession in accessions])
    with open(projects_file, "a+") as pfh:
        pfh.writelines([bioproject + "\n" for bioproject in new_projects])

    # TODO: introduce page numbering for when list gets long
