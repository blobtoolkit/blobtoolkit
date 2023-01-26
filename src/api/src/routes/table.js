const Field = require("../models/field");
const utils = require("../shared/functions/utils");

const expandKeys = (value, keys, slot) => {
  let hasCategory = keys && keys.length > 0;
  let isArray = Array.isArray(value);
  if (isArray && hasCategory && !isNaN(slot)) {
    for (let v of value) {
      v[slot] = keys[v[slot]];
    }
    return value;
  }
  if (hasCategory && !isArray) {
    return keys[value];
  }
  return value;
};

const getTable = async (req, res) => {
  let fields = req.params.field_ids.split(",");
  if (!fields.includes("identifiers")) {
    fields = ["identifiers", ...fields];
  }
  let table = [fields];
  let success = true;
  let windowSize = req.query.window_size;
  for (let field_id of fields) {
    let data;
    if (windowSize) {
      let field = new Field(
        `${field_id}_windows${windowSize == 0.1 ? "" : `_${windowSize}`}`,
        {
          id: req.params.dataset_id,
        }
      );
      data = await field.loadData();
    }
    if (!data) {
      let field = new Field(field_id, { id: req.params.dataset_id });
      data = await field.loadData();
    }
    console.log({ field_id, data });
    if (data && data.hasOwnProperty("values")) {
      if (table.length == 1) {
        for (let id of data.values) {
          table.push([id]);
        }
      } else if (data.values.length == table.length - 1) {
        for (let i = 0; i < data.values.length; i++) {
          table[i + 1].push(
            expandKeys(data.values[i], data.keys, data.category_slot)
          );
        }
      } else {
        success = false;
        break;
      }
    } else {
      success = false;
      break;
    }
  }
  if (success) {
    res.setHeader("content-type", "application/json");
    res.json(table);
  } else {
    res.sendStatus(404);
  }
};

const getTableByIndex = async (req, res) => {
  res.setHeader("content-type", "application/json");
  let field = new Field(req.params.field_id, { id: req.params.dataset_id });
  let data = await field.loadDataAtIndex(req.params.index);
  if (data) {
    res.json(data);
  } else {
    res.sendStatus(404);
  }
};

module.exports = {
  getTable,
  getTableByIndex,
};
