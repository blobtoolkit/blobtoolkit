const config = require("../../src/config/main");
const checkCompression = require("../functions/io").checkCompression;

const getImage = async (req, res) => {
  file = await checkCompression(
    config.filePath + "/" + req.params.dataset_id + "/blob.circle.png"
  );
  if (file) {
    res.setHeader("content-type", "image/png");
    if (file.endsWith(".gz")) {
      res.setHeader("content-encoding", "gzip");
    }
    res.sendFile(file);
  } else {
    res.sendStatus(404);
  }
};
const getImageView = async (req, res) => {
  let type = "";
  if (req.params.view == "blob") {
    type = ".circle";
  }
  if (req.query.format == "svg") {
    let file =
      config.filePath +
      "/" +
      req.params.dataset_id +
      "/" +
      req.params.view +
      type +
      ".svg";
    file = await checkCompression(file);
    if (file) {
      res.setHeader("content-type", "image/svg+xml");
      if (file.endsWith(".gz")) {
        res.setHeader("content-encoding", "gzip");
      }
      res.sendFile(file);
    } else {
      res.sendStatus(404);
    }
  } else {
    res.setHeader("content-type", "image/png");
    let file =
      config.filePath +
      "/" +
      req.params.dataset_id +
      "/" +
      req.params.view +
      type +
      ".png";
    file = await checkCompression(file);
    if (file) {
      res.setHeader("content-type", "image/png");
      if (file.endsWith(".gz")) {
        res.setHeader("content-encoding", "gzip");
      }
      res.sendFile(file);
    } else {
      res.sendStatus(404);
    }
  }
};
const getImageType = async (req, res) => {
  let type = req.params.type;
  if (type == "none") type = "circle";
  if (req.query.format == "svg") {
    let file =
      config.filePath + "/" + req.params.dataset_id + "/blob." + type + ".svg";
    file = await checkCompression(file);
    if (file) {
      res.setHeader("content-type", "image/svg+xml");
      if (file.endsWith(".gz")) {
        res.setHeader("content-encoding", "gzip");
      }
      res.sendFile(file);
    } else {
      res.sendStatus(404);
    }
  } else {
    let file =
      config.filePath + "/" + req.params.dataset_id + "/blob." + type + ".png";
    file = await checkCompression(file);
    if (file) {
      res.setHeader("content-type", "image/png");
      if (file.endsWith(".gz")) {
        res.setHeader("content-encoding", "gzip");
      }
      res.sendFile(file);
    } else {
      res.sendStatus(404);
    }
  }
};

module.exports = { getImage, getImageView, getImageType };
