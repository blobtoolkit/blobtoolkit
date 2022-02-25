import { createAction, handleAction, handleActions } from "redux-actions";
import { getAvailableDatasetIds, getAvailableDatasets } from "./repository";

import { createSelector } from "reselect";
import { format as d3Format } from "d3-format";
import store from "../store";

export const setDatasetPage = createAction("SET_DATASET_PAGE");

export const datasetPage = handleAction(
  "SET_DATASET_PAGE",
  (state, action) => action.payload,
  0
);

export const getDatasetPage = (state) => state.datasetPage;

export const setDatasetPageSize = createAction("SET_DATASET_PAGESIZE");

export const datasetPageSize = handleAction(
  "SET_DATASET_PAGESIZE",
  (state, action) => action.payload,
  10
);

export const getDatasetPageSize = (state) => state.datasetPageSize;

export const setDatasetSorted = createAction("SET_DATASET_SORTED");

export const datasetSorted = handleAction(
  "SET_DATASET_SORTED",
  (state, action) => action.payload,
  []
);

export const getDatasetSorted = (state) => state.datasetSorted;

export const setDatasetColumns = createAction("SET_DATASET_COLUMNS");

export const datasetColumns = handleAction(
  "SET_DATASET_COLUMNS",
  (state, action) => action.payload,
  [
    "id",
    "taxon",
    "accession",
    "-hits-count",
    "-hits-span",
    "-hits-n50",
    "-busco-lineage",
    "-busco-string",
    "reads",
  ]
);

export const getDatasetColumns = (state) => state.datasetColumns;

const requireTotal = (d, key) => {
  if (d.summaryStats && d.summaryStats.hits && d.summaryStats.hits.total) {
    return d.summaryStats.hits.total[key];
  }
  return "";
};

export const buscoFields = {
  c: "Complete",
  d: "Duplicated",
  f: "Fragmented",
  s: "Single-copy",
  m: "Missing",
  n: "Total",
  string: "String",
  lineage: "Lineage",
};

const requireBusco = (d, lineage, field) => {
  if (d.summaryStats && d.summaryStats.busco) {
    if (!lineage && Object.keys(d.summaryStats.busco).length > 0) {
      lineage = Object.keys(d.summaryStats.busco)[0];
    }
    if (lineage) {
      if (field == "lineage") {
        return lineage;
      }
      if (d.summaryStats.busco[lineage]) {
        return d.summaryStats.busco[lineage][field];
      }
    }
    return "";
  }
  return "";
};

const requireComposition = (d, key) => {
  let pct = (n) => d3Format(".1f")(n * 100);
  if (
    d.summaryStats &&
    d.summaryStats.baseComposition &&
    d.summaryStats.baseComposition[key]
  ) {
    return pct(d.summaryStats.baseComposition[key]);
  }
  return "";
};

const requireHits = (d, taxon, key) => {
  if (!taxon) taxon = "total";
  if (
    d.summaryStats &&
    d.summaryStats.hits &&
    d.summaryStats.hits[taxon] &&
    d.summaryStats.hits[taxon][key]
  ) {
    return d.summaryStats.hits[taxon][key];
  }
  return "";
};

const requireMapping = (d) => {
  if (d.summaryStats && d.summaryStats.readMapping) {
    return Object.keys(d.summaryStats.readMapping).length;
  }
  return "";
};

const requireTaxonomy = (d, key) => {
  if (
    d.summaryStats &&
    d.summaryStats.taxonomy &&
    d.summaryStats.taxonomy[key]
  ) {
    return d.summaryStats.taxonomy[key];
  }
  return "";
};

const requireStats = (d, key) => {
  let pct = (n) => d3Format(".1f")(n * 100);
  let rat = d3Format(".2r");
  if (d.summaryStats && d.summaryStats.stats) {
    if (key == "n50overSpan") {
      if (d.summaryStats.stats[key]) {
        return rat(d.summaryStats.stats[key]);
      } else if (d.summaryStats.stats.spanOverN50) {
        return rat(1 / d.summaryStats.stats.spanOverN50);
      }
      return "";
    } else if (d.summaryStats.stats[key]) {
      return pct(d.summaryStats.stats[key]);
    }
    return "";
  }
  return "";
};

let colInfo = {
  id: {
    title: "Dataset ID",
    func: (d) => d.id,
    group: "Identifiers",
  },
  taxon: {
    title: "Taxon",
    func: (d) => d.taxon_name,
    width: 300,
    group: "Identifiers",
  },
  accession: {
    title: "Accession",
    func: (d) => d.accession,
    width: 150,
    group: "Identifiers",
  },
  count: {
    title: "Sequences",
    func: (d, taxon) => requireHits(d, taxon, "count"),
    numeric: true,
    group: "Assembly statistics",
  },
  span: {
    title: "Span (bp)",
    func: (d, taxon) => requireHits(d, taxon, "span"),
    numeric: true,
    group: "Assembly statistics",
  },
  n50: {
    title: "N50 (bp)",
    func: (d, taxon) => requireHits(d, taxon, "n50"),
    numeric: true,
    group: "Assembly statistics",
  },
  busco: {
    title: (d, lineage, field) => "Closest " + field,
    func: (d, lineage, field) => requireBusco(d, lineage, field),
    group: "BUSCO",
  },
  gc: {
    title: "GC (%)",
    func: (d) => requireComposition(d, "gc"),
    numeric: true,
    group: "Assembly statistics",
  },
  at: {
    title: "AT (%)",
    func: (d) => requireComposition(d, "at"),
    numeric: true,
    group: "Assembly statistics",
  },
  n: {
    title: "N (%)",
    func: (d) => requireComposition(d, "n"),
    numeric: true,
    group: "Assembly statistics",
  },
  lineage: {
    title: "Lineage",
    func: (d) => requireTaxonomy(d, "lineage"),
    width: 500,
    group: "Taxonomy",
  },
  taxid: {
    title: "Taxon ID",
    func: (d) => requireTaxonomy(d, "taxid"),
    group: "Taxonomy",
  },
  rank: {
    title: "Target Rank",
    func: (d) => requireTaxonomy(d, "targetRank"),
    group: "Taxonomy",
  },
  target: {
    title: "Target taxon",
    func: (d) => requireTaxonomy(d, "target"),
    width: 150,
    group: "Taxonomy",
  },
  nohit: {
    title: "No hit (%)",
    func: (d) => requireStats(d, "noHit"),
    numeric: true,
    group: "Summary statistics",
  },
  hit: {
    title: "Hits matching target (%)",
    func: (d) => requireStats(d, "target"),
    numeric: true,
    group: "Summary statistics",
  },
  ratio: {
    title: "N50:Span ratio",
    func: (d) => requireStats(d, "n50overSpan"),
    numeric: true,
    group: "Summary statistics",
  },
  reads: {
    title: "Read sets",
    func: (d) => requireMapping(d),
    numeric: true,
    group: "Coverage",
  },
};

export const datasetSummaries = createSelector(
  getAvailableDatasetIds,
  getAvailableDatasets,
  getDatasetColumns,
  (ids, datasets, colNames) => {
    let data = [];
    ids.forEach((id) => {
      let values = {};
      colNames.forEach((key) => {
        if (key.match("-busco-")) {
          let [lineage, field] = key.split("-busco-");
          colInfo[key] = Object.assign({}, colInfo["busco"]);
          if (lineage) {
            if (!lineage.endsWith("_odb9")) lineage += "_odb9";
            colInfo[key].title = lineage + " " + buscoFields[field];
          } else {
            colInfo[key].title = colInfo["busco"].title(
              datasets[id],
              lineage,
              field
            );
          }
          values[key] = colInfo["busco"].func(datasets[id], lineage, field);
          if (field == "string") {
            colInfo[key].width = 300;
          } else if (field == "lineage") {
            colInfo[key].width = 150;
          }
        } else if (key.match("-hits-")) {
          let [taxon, field] = key.split("-hits-");
          colInfo[key] = Object.assign({}, colInfo[field]);
          if (taxon)
            colInfo[key].title = colInfo[field].title + " (" + taxon + ")";
          values[key] = colInfo[field].func(datasets[id], taxon, field);
        } else {
          values[key] = colInfo[key].func(datasets[id]);
        }
      });
      data.push(values);
    });
    return data;
  }
);

export const hitTaxa = createSelector(
  getAvailableDatasetIds,
  getAvailableDatasets,
  (ids, datasets) => {
    let taxa = {};
    ids.forEach((id) => {
      if (datasets[id].summaryStats && datasets[id].summaryStats.hits) {
        Object.keys(datasets[id].summaryStats.hits).forEach((taxon) => {
          if (taxon != "total") {
            taxa[taxon] = true;
          }
        });
      }
    });
    return Object.keys(taxa).sort((a, b) => (a > b ? 1 : -1));
  }
);

export const buscoLineages = createSelector(
  getAvailableDatasetIds,
  getAvailableDatasets,
  (ids, datasets) => {
    let lineages = {};
    ids.forEach((id) => {
      if (datasets[id].summaryStats && datasets[id].summaryStats.busco) {
        Object.keys(datasets[id].summaryStats.busco).forEach((lineage) => {
          lineages[lineage] = true;
        });
      }
    });
    return Object.keys(lineages).sort((a, b) => (a > b ? 1 : -1));
  }
);

const ascComma = (a, b) => {
  if (a.length === b.length) {
    return a > b ? 1 : -1;
  }
  return a.length > b.length ? 1 : -1;
};

export const columnInfo = createSelector(
  () => colInfo,
  getDatasetColumns,
  datasetSummaries,
  (obj, columns) => {
    let info = {};
    Object.keys(obj).forEach((key) => {
      let group = obj[key].group;
      if (!info[group]) {
        info[group] = {};
      }
      if (key.startsWith("-hits-")) {
        let base = key.replace("-hits-", "");
        delete info[group][base];
      }
      info[group][key] = {
        title: obj[key].title,
      };
      if (columns.includes(key)) {
        info[group][key].active = true;
      }
    });
    return info;
  }
);

export const listingColumns = createSelector(columnInfo, (groupedColumns) => {
  let columns = [];
  Object.keys(groupedColumns).forEach((groupName) => {
    let group = groupedColumns[groupName];
    let cols = [];
    Object.keys(group).forEach((key) => {
      if (group[key].active) {
        let config = {};
        config.id = key;
        config.Header = colInfo[key].title;
        config.accessor = (d) => {
          if (colInfo[key].numeric) {
            return d[key].toLocaleString();
          }
          return d[key];
        };
        config.width = colInfo[key].width ? colInfo[key].width : 100;
        if (colInfo[key].numeric) {
          config.style = { textAlign: "right" };
          config.sortMethod = ascComma;
        }
        cols.push(config);
      }
    });
    if (cols.length > 0) {
      let config = {};
      config.Header = groupName;
      config.columns = cols;
      columns.push(config);
    }
  });
  // colNames.forEach(key=>{
  //   let config = {}
  //   config.id = key
  //   config.Header = colInfo[key].title
  //   config.accessor = d => {
  //     if (colInfo[key].numeric){
  //       return d[key].toLocaleString()
  //     }
  //     return d[key]
  //   }
  //   config.width = colInfo[key].width ? colInfo[key].width : 100
  //   if (colInfo[key].numeric){
  //     config.style = {textAlign:'right'}
  //     config.sortMethod = ascComma
  //   }
  //   columns.push(config)
  // })
  return columns;
});

export const getDatasetCSVdata = createSelector(
  columnInfo,
  datasetSummaries,
  (columns, data) => {
    let arr = [];
    let cells = [];
    Object.keys(columns).forEach((groupName) => {
      let group = columns[groupName];
      Object.keys(group).forEach((key) => {
        if (group[key].active) {
          cells.push(colInfo[key].title);
        }
      });
    });
    arr.push(cells);
    data.forEach((row) => {
      cells = [];
      Object.keys(columns).forEach((groupName) => {
        let group = columns[groupName];
        Object.keys(group).forEach((key) => {
          if (group[key].active) {
            cells.push(row[key]);
          }
        });
      });
      arr.push(cells);
    });
    return arr;
  }
);

export const datasetTableReducers = {
  datasetPage,
  datasetPageSize,
  datasetSorted,
  datasetColumns,
};
