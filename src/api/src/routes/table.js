const Field = require("../models/field");
const utils = require("../shared/functions/utils");

const expandKeys = (value, keys, slot) => {
  let hasCategory = keys && keys.length > 0;
  let isArray = Array.isArray(value);
  if (isArray && hasCategory && !isNaN(slot) && !isNaN(slot)) {
    if (Array.isArray(value[0])) {
      for (let v of value) {
        v[slot] = keys[v[slot]];
      }
    } else {
      value[slot] = keys[value[slot]];
    }
    return value;
  }
  if (hasCategory && !isArray) {
    return keys[value];
  }
  return value;
};

const add_window_to_table = ({
  table,
  indices,
  headers,
  field_id,
  data,
  hasWindows,
}) => {
  data.values = data.values.slice(5, 10);
  if (table.length == 1) {
    for (let i = 0; i < data.values.length; i++) {
      indices[i] = [table.length];
      indices.state = "default";
      table.push([data.values[i]]);
    }
    return { success: true, table };
  }
  if (hasWindows) {
    let prevTable;
    let start = 0;
    if (indices.state == "default") {
      prevTable = JSON.parse(JSON.stringify(table));
      table = [table[0]];
      headers.splice(1, 1, "start", "end");
    } else if (data.headers.length && data.headers.length > 1) {
      let index = headers.indexOf(field_id);
      let new_headers = data.headers.map((h) =>
        h.startsWith(field_id) ? h : `${field_id}_${h}`
      );
      headers.splice(index, 1, ...new_headers);
    }
    let j = 0;
    for (let i = 0; i < data.values.length; i++) {
      let rowIndex = indices[i][0];
      for (let value of data.values[i]) {
        j += 1;
        if (indices.state == "default") {
          table[j] = [...prevTable[rowIndex]];
          table[j].push(start);
          start += value[0];
          table[j].push(start);
        } else {
          if (Array.isArray(value)) {
            let v = expandKeys(value, data.keys, data.category_slot);
            table[j].push(...v);
          } else {
            table[j].push(value);
          }
        }
      }
    }
    indices.state = "expanded";
    return { success: true, table };
  }
  for (let i = 0; i <= data.values.length; i++) {
    for (let index of indices[i]) {
      table[index].push([data.values[i]]);
    }
  }
  return { success: true, table };
};

const add_to_table = (table, data) => {
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
    return false;
  }
  return true;
};

const getTable = async (req, res) => {
  let fields = req.params.field_ids.split(",");
  let windowSize = req.query.window_size;
  let indices;
  if (windowSize) {
    fields = ["length", ...fields];
    indices = [];
  }
  if (!fields.includes("identifiers")) {
    fields = ["identifiers", ...fields];
  }
  let table = [fields];
  let headers = [...fields];
  let success = true;
  let hasWindows = false;

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
    if (data) {
      hasWindows = true;
    } else {
      let field = new Field(field_id, { id: req.params.dataset_id });
      data = await field.loadData();
    }
    if (data && data.hasOwnProperty("values")) {
      if (windowSize) {
        ({ success, table } = add_window_to_table({
          table,
          indices,
          headers,
          field_id,
          data,
          hasWindows,
        }));
      } else {
        success = add_to_table(table, data);
      }
    } else {
      success = false;
      break;
    }
  }
  if (success) {
    table[0] = headers;
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
