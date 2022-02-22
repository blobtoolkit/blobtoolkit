const Repository = require("../shared/models/repository");
const config_main = require("../../src/config/main");
const io = require("../functions/io");

module.exports = Repository;

Repository.prototype.loadMeta = () => {
  let path = config_main.filePath;
  return io.readJSON(path + "/meta.json");
};
