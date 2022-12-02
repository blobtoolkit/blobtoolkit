const checkCompression = require("../functions/io").checkCompression;

const config = require("../config/main");
const getSummaryById = async (req, res) => {
  file = await checkCompression(
    config.filePath + "/" + req.params.dataset_id + "/summary.json"
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
  getSummaryById,
};
