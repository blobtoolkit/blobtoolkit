#!/usr/bin/env node

import { csvToObj } from "csv-to-js-parser";
import fs from "fs";

// Generate a targets.json file for the BlobToolKit Viewer target tree
// TODO: replace hardcoded file paths with fetch commands

// TODO: replace idFilePath with fetch from BlobToolKit
const idFilePath =
  "~/projects/blobtoolkit/blobtoolkit/src/api/scripts/validated_accession_id.tsv";

// TODO: replace csvFilePath with fetch from https://goat.genomehubs.org/search?taxonomy=ncbi&query=tax_tree%282759%5BEukaryota%5D%29&result=assembly&includeEstimates=false&summaryValues=count&offset=0&fields=assembly_type&names=&ranks=subspecies%2Cspecies%2Cgenus%2Cfamily%2Corder%2Cclass%2Cphylum%2Ckingdom%2Csuperkingdom&size=100#tax_tree(2759%5BEukaryota%5D)
const csvFilePath =
  "~/projects/blobtoolkit/blobtoolkit/src/api/scripts/goat_assembly.tsv";

const outFilePath =
  "~/projects/blobtoolkit/blobtoolkit/src/api/scripts/targets.json";

const ranks = [
  "superkingdom",
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species",
];

let ids = {};

let parents = {};

let nodes = {};

fs.readFileSync(idFilePath, "utf-8")
  .split(/\r?\n/)
  .forEach((line) => {
    let [accession, id] = line.replaceAll(/"/g, "").split("\t");
    ids[accession] = id;
    accession = accession.replace(/\.\d+/, "");
    ids[accession] = id;
  });

const data = fs.readFileSync(csvFilePath).toString();

let rows = csvToObj(data, "\t");

let n = 0;
nodes["root"] = { n, r: "root", ta: 0, ts: 0, d: {} };
n++;
for (let row of rows) {
  let speciesLineage = ranks.map((rank) => row[rank]).join("");
  let newSpecies = !nodes[speciesLineage];
  let parent = "root";
  let lineage = "";
  for (let rank of ranks) {
    let prevLineage = lineage;
    lineage += row[rank];
    if (!nodes[lineage]) {
      nodes[lineage] = { n, r: row[rank], ta: 0, ts: 0, d: {} };
      parents[lineage] = {
        parent,
        name:
          row[rank] ||
          `Other ${parents[prevLineage].name}`.replace("Other Other", "Other"),
      };
      n++;
    }
    parent = lineage;
    if (newSpecies) {
      nodes[lineage].ts += 1;
    }
    nodes[lineage].ta += 1;
    if (rank == "species") {
      let id =
        ids[row.assembly_id] || ids[row.assembly_id.replace(/\.\d+/, "")];
      nodes[lineage].d[id] = { n };
      n++;
    }
  }
}

for (let [lineage, obj] of Object.entries(parents)) {
  nodes[obj.parent].d[obj.name] = nodes[lineage];
  nodes[obj.parent].d[obj.name] = nodes[lineage];
}

for (let obj of Object.values(nodes.root.d)) {
  nodes.root.ta += obj.ta;
  nodes.root.ts += obj.ts;
}

fs.writeFileSync(outFilePath, JSON.stringify(nodes.root));
