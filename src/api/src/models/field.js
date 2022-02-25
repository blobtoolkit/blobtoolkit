const config = require("../../src/config/main");
const io = require("../functions/io");
const Filter = require("../shared/models/filter");
const Field = require("../shared/models/field");

module.exports = Field;

const loadData = async function (useAPI) {
  if (this.data) return Promise.resolve(this.data);
  let filePath = this.filePath || config.filePath;
  this.data = await io.readJSON(
    filePath + "/" + this.dataset.id + "/" + this._id + ".json"
  );
  return Promise.resolve(this.data);
};

const loadDataAtIndex = async function (index) {
  let data = await this.loadData();
  if (!data) return undefined;
  let lastIndex = data.values.length - 1;
  let values = [];
  let keys = false;
  let ret;
  if (data.hasOwnProperty("keys")) {
    keys = data.keys;
  }
  String(index)
    .split(",")
    .forEach((r) => {
      let indices = r.split("-");
      if (r.match("-") && (!indices[0] || !indices[1])) {
        ret = Promise.resolve(undefined);
      }
      indices[0] = indices[0] ? indices[0] * 1 : 0;
      indices[1] = indices[1] ? indices[1] * 1 : indices[0];
      if (indices[1] < indices[0] || indices[1] > lastIndex) {
        ret = Promise.resolve(undefined);
      } else {
        for (var i = indices[0]; i <= indices[1]; i++) {
          if (keys && keys.length > 0) {
            if (typeof data.values[i] == "object") {
              values.push(data.values[i]);
            } else {
              values.push(data.keys[data.values[i]]);
            }
          } else {
            values.push(data.values[i]);
          }
        }
      }
    });
  if (keys && keys.length > 0) {
    ret = ret || Promise.resolve({ values, keys });
  } else {
    ret = ret || Promise.resolve(values);
  }
  return ret;
};

Field.prototype.loadData = loadData;
Field.prototype.loadDataAtIndex = loadDataAtIndex;
