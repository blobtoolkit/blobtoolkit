const appRoot = require("app-root-path");
const Dataset = require("../models/dataset");
const Field = require("../models/field");
const utils = require("../shared/functions/utils");

module.exports = function (app, db) {
  var directory = app.locals.directory;

  function values(ret, meta, path, index) {
    meta.forEach(function (n) {
      var obj = {};
      if (n.children) {
        obj[n.id] = [];
        ret.push(obj);
        values(obj[n.id], n.children, path, index);
      } else {
        var js = require(path + n.id + ".json");
        obj[n.id] = js[index];
        ret.push(obj);
      }
    });
  }

  function slice(index, types, meta, path) {
    var obj = { index: index };
    types.forEach(function (t) {
      obj[t] = [];
      values(obj[t], meta[t], path, index);
    });
    return obj;
  }

  function list_fields(arr, dataset_id, index, ret) {
    arr.forEach(async (field_meta) => {
      let field = new Field(field_meta.id, { id: dataset_id });
      if (field) {
        let data = await field.loadDataAtIndex(index);
        if (data) {
          ret[field_meta.id] = data;
        }
        if (field_meta.children) {
          ret = Object.assign(
            ret,
            list_fields(field_meta.children, dataset_id, index, ret)
          );
        }
        if (field_meta.data) {
          ret = Object.assign(
            ret,
            list_fields(field_meta.data, dataset_id, index, ret)
          );
        }
      }
    });
    return ret;
  }

  app.get("/api/v1/:dataset_id/slice/:index", async (req, res) => {
    // index can be single value, comma separated, or range
    // open-ended ranges will be expanded to end of dataset
    let dataset = new Dataset(req.params.dataset_id);
    let meta = await dataset.loadMeta();
    let ret = list_fields(
      meta.fields,
      req.params.dataset_id,
      req.params.index,
      {}
    );

    console.log(ret);
    res.json(ret);
  });
};
