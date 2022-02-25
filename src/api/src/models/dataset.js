const config = require("../../src/config/main");
const io = require("../functions/io");
const Field = require("./field");
const utils = require("../shared/functions/utils");
const Dataset = require("../shared/models/dataset");

module.exports = Dataset;

const loadBlobDB = async function (file = this.blobDBFile) {
  if (this.blobDB) return Promise.resolve(this.blobDB);
  if (!file)
    return Promise.reject(Error("Cannot loadBlobDB without a filename"));
  this.blobDB = await io.readJSON(file);
  this.blobDBFile = file;
  return Promise.resolve(this.blobDB);
};

const prepareMeta = async function () {
  if (this.meta) return Promise.resolve(this.meta);
  let json = await this.loadBlobDB();
  let meta = {};
  meta.filePath = config.filePath;
  meta.outFilePath = config.outFilePath;
  meta.id = json.title
    .replace(/blobdb/i, "")
    .replace(/json/i, "")
    .replace(/[\W]+/, "");
  meta.name = json.title;
  meta.description = json.description || "imported from " + meta.name;
  meta.records = json.seqs;
  meta.record_type = "contigs";
  let template = {};
  let values = {};
  let fields = [];
  meta.fields = fields;
  fields.push(
    this.addField("gc", {
      name: "GC",
      description: "per contig GC percentage",
      type: "variable",
      datatype: "float",
      range: [0, 1],
      scale: "scaleLinear",
      preload: true,
      blobDBpath: ["gc"],
    })
  );
  fields.push(
    this.addField("length", {
      name: "Length",
      description: "per contig length",
      type: "variable",
      datatype: "integer",
      range: [1, json.length],
      scale: "scaleLog",
      preload: true,
      blobDBpath: ["length"],
    })
  );
  if (json.n_count > 0) {
    fields.push(
      this.addField("ncount", {
        name: "N-count",
        description: "Ns per contig",
        type: "variable",
        datatype: "integer",
        range: [1, json.n_count],
        scale: "scaleLog",
        blobDBpath: ["n_count"],
      })
    );
  }
  let covLibs = [];
  let readCovLibs = [];
  Object.keys(json.covLibs).forEach((key) => {
    let lib = json.covLibs[key];
    let info = { name: lib.name, blobDBpath: ["covs", key] };
    if (key == "cov0") {
      info["preload"] = true;
    }
    covLibs.push(this.addField(key + "_cov", info));
    readCovLibs.push(
      this.addField(key + "_read_cov", {
        name: lib.name,
        blobDBpath: ["read_cov", key],
      })
    );
  });
  fields.push(
    this.addField("covs", {
      name: "Coverage",
      description: "coverage per contig",
      type: "variable",
      datatype: "float",
      range: [1, 100000],
      scale: "scaleLog",
      children: covLibs,
    })
  );
  fields.push(
    this.addField("read_cov", {
      name: "Read coverage",
      description: "read coverage per contig",
      type: "variable",
      datatype: "float",
      range: [1, 100000],
      scale: "scaleLog",
      children: readCovLibs,
    })
  );
  let hitLibs = [];
  Object.keys(json.hitLibs).forEach((key) => {
    let lib = json.hitLibs[key];
    hitLibs.push(
      this.addField(key, { name: lib.name, blobDBpath: ["hits", key] })
    );
  });
  fields.push(
    this.addField("hits", {
      name: "Hits",
      description: "Blast hit taxonomy IDs",
      type: "hit",
      datatype: "array",
      lookup: "ncbi_taxonomy",
      children: hitLibs,
    })
  );
  let taxrules = [];
  let levels = [
    "superkingdom",
    "phylum",
    "order",
    "family",
    "genus",
    "species",
  ];
  let taxfields = [];
  json.taxrules.forEach((rule) => {
    let taxlevels = [];
    levels.forEach((level) => {
      let data = [];
      data.push(
        this.addField(rule + "_" + level + "_score", {
          name: "score",
          type: "variable",
          datatype: "float",
          range: [1, 10000],
          scale: "scaleSqrt",
          blobDBpath: ["taxonomy", rule, level, "score"],
          preload: false,
          active: false,
        })
      );
      data.push(
        this.addField(rule + "_" + level + "_cindex", {
          name: "c_index",
          type: "variable",
          datatype: "integer",
          range: [1, 10],
          scale: "scaleLinear",
          blobDBpath: ["taxonomy", rule, level, "c_index"],
          preload: false,
          active: false,
        })
      );
      let info = {
        name: level,
        data: data,
        blobDBpath: ["taxonomy", rule, level, "tax"],
      };
      if (rule + "_" + level == "bestsumorder_phylum") {
        info["preload"] = true;
      }
      taxlevels.push(this.addField(rule + "_" + level, info));
    });
    taxrules.push(this.addField(rule, { children: taxlevels }));
  });
  fields.push(
    this.addField("taxonomy", {
      name: "Taxonomy",
      description: "BLAST-assigned taxonomy",
      type: "category",
      datatype: "string",
      children: taxrules,
    })
  );

  this.meta = meta; //await waitOn(meta,json);
  return Promise.resolve(this.meta);
};

const storeMeta = async function () {
  if (!this.meta)
    return Promise.reject(Error("Cannot storeMeta if meta is undefined"));
  let filePath = this.outFilePath || config.outFilePath;
  utils.removeNestedKeys(
    this.meta.fields,
    ["dataset", "filters", "scale"],
    ["_data", "_children"]
  );
  this.meta.fields = this.meta.fields.map((field) => {
    let obj = {};
    Object.keys(field).forEach((key) => {
      obj[key] = field[key];
    });
    return obj;
  });
  let success = await io.writeJSON(
    filePath + "/" + this.id + "/meta.json",
    this.meta
  );
  return success;
};

const defaultLinks = {
  links: {
    taxon: {
      taxid: {
        ENA: "https://www.ebi.ac.uk/ena/browser/view/Taxon:{taxid}",
      },
      species: {
        Wikipedia: "https://wikipedia.org/wiki/{species}",
      },
    },
    blobtoolkit: {
      commit: {
        Github: "{pipeline}/tree/{commit}",
      },
    },
    assembly: {
      biosample: {
        ENA: "https://www.ebi.ac.uk/ena/browser/view/{biosample}",
      },
      bioproject: {
        ENA: "https://www.ebi.ac.uk/ena/browser/view/{bioproject}",
      },
    },
    position: [
      {
        patterns: [
          {
            title: "NCBI RefSeq",
            template: "https://www.ncbi.nlm.nih.gov/nuccore/{subject}",
            regex: "^([NXW][A-Z]_.+.d+)$",
          },
          {
            title: "UniProt",
            template: "https://www.uniprot.org/uniprot/{subject}",
            regex:
              "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})",
          },
          {
            title: "ENA",
            template: "https://www.ebi.ac.uk/ena/browser/view/{subject}",
            regex: "^(.+)$",
          },
        ],
      },
      {
        UniProt: "https://www.uniprot.org/uniprot/{subject}",
      },
    ],
    position_reverse: "_aa_",
  },
};

const loadMeta = async function () {
  if (!this.meta) {
    let filePath = this.filePath || config.filePath;
    let defaultMeta = await io.readJSON(filePath + "/default.json");
    if (!defaultMeta && config.use_default) defaultMeta = defaultLinks;
    this.meta = await io.readJSON(filePath + "/" + this.id + "/meta.json");
    if (!this.meta) this.meta = {};
    if (
      this.meta.record_type == "record" &&
      this.meta.assembly &&
      this.meta.assembly.level
    ) {
      this.meta.record_type = this.meta.assembly.level;
    }
    if (!this.meta.hasOwnProperty("revision")) {
      this.meta.revision = 0;
    }
    if (defaultMeta) {
      Object.keys(defaultMeta).forEach((key) => {
        if (!this.meta.hasOwnProperty(key)) {
          this.meta[key] = defaultMeta[key];
        } else {
          Object.keys(defaultMeta[key]).forEach((k) => {
            if (!this.meta[key].hasOwnProperty(k)) {
              this.meta[key][k] = defaultMeta[key][k];
            }
          });
        }
      });
    }
  }
  if (this.meta.hasOwnProperty("id")) {
    this.addFields(this.meta.fields);
  }

  return Promise.resolve(this.meta);
};

const storeLineages = async function () {
  let lineages = this.blobDB.lineages;
  let levels = Object.keys(lineages[Object.keys(lineages)[0]]);
  let values = {};
  Object.keys(lineages).forEach((taxid) => {
    (values[taxid] = []),
      levels.forEach((level, i) => {
        values[taxid][i] = lineages[taxid][level] || undefined;
      });
  });
  let filePath = this.outFilePath || config.outFilePath;
  let success = await io.writeJSON(
    filePath + "/" + this.id + "/lineages.json",
    { keys: levels, values: values }
  );
  return success;
};

const storeValues = async function (id, field) {
  field =
    field ||
    this.meta.fields.find(function (f) {
      return f._id == id;
    });
  let values = [];
  let successes = [];
  if (field._children) {
    field._children.forEach(async (child) => {
      let success = this.storeValues(child._id, child);
      successes.push(success);
    });
  } else {
    this.blobDB.order_of_blobs.forEach(async (contig) => {
      let val = utils.valueAtPath(
        this.blobDB.dict_of_blobs[contig],
        field._blobDBpath
      );
      if (field._blobDBpath[0] == "hits") {
        value = await utils.groupValuesBy(val || [], "taxId", "score");
      } else {
        value = val;
      }
      values.push(value);
    });
    let filePath = this.outFilePath || config.outFilePath;
    if (
      field._blobDBpath[0] == "taxonomy" &&
      field._blobDBpath[field._blobDBpath.length - 1] == "tax"
    ) {
      let success;
      let result = utils.asKeyValue(values);
      success = io.writeJSON(
        filePath + "/" + this.id + "/" + field._id + ".json",
        result
      );
      successes.push(success);
    } else {
      let result = { values: values };
      let success = io.writeJSON(
        filePath + "/" + this.id + "/" + field._id + ".json",
        result
      );
      successes.push(success);
    }
  }
  if (field._data) {
    field._data.forEach((child) => {
      let success = this.storeValues(child._id, child);
      successes.push(success);
    });
  }
  return Promise.all(successes);
};

const storeAllValues = function () {
  let successes = [];
  this.meta.fields.forEach(async (field) => {
    let out = this.storeValues(field._id);
    successes = successes.concat(out);
  });
  let filePath = this.outFilePath || config.outFilePath;
  let success = io.writeJSON(
    filePath + "/" + this.id + "/identifiers.json",
    this.blobDB.order_of_blobs
  );
  successes.push(success);
  return Promise.all(successes);
};

Dataset.prototype.loadBlobDB = loadBlobDB;
Dataset.prototype.prepareMeta = prepareMeta;
Dataset.prototype.loadMeta = loadMeta;
Dataset.prototype.storeMeta = storeMeta;
Dataset.prototype.storeLineages = storeLineages;
Dataset.prototype.storeValues = storeValues;
Dataset.prototype.storeAllValues = storeAllValues;
