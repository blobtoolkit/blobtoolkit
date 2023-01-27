const Field = require("../models/field");
const utils = require("../shared/functions/utils");

const getFieldById = async (req, res) => {
  res.setHeader("content-type", "application/json");
  //res.setHeader('Access-Control-Allow-Origin', 'http://localhost:8080');
  let windowSize = req.query.window_size;
  if (windowSize) {
    let field = new Field(
      `${req.params.field_id}_windows${
        windowSize == 0.1 ? "" : `_${windowSize}`
      }`,
      {
        id: req.params.dataset_id,
      }
    );
    data = await field.loadData();
  }
  if (data) {
    hasWindows = true;
  } else {
    let field = new Field(req.params.field_id, { id: req.params.dataset_id });
    data = await field.loadData();
  }
  if (data && data.hasOwnProperty("values")) {
    res.json(data);
  } else {
    res.sendStatus(404);
  }
};

const getFieldByIdIndex = async (req, res) => {
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
  getFieldById,
  getFieldByIdIndex,
};
