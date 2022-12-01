const checkCompression = require("../functions/io").checkCompression;
const config = require("../../src/config/main");

const getIdentifiers = async (req, res) => {
  file = await checkCompression(
    config.filePath + "/" + req.params.dataset_id + "/identifiers.json"
  );
  if (file) {
    res.setHeader("content-type", "application/json");
    if (file.endsWith(".gz")) {
      res.setHeader("content-encoding", "gzip");
    }
    res.sendFile(file);
  } else {
    res.sendStatus(404);
  }
};

module.exports = {
  getIdentifiers,
};
