const Promise = require("promise");
const fs = require("fs");
const zlib = require("zlib");
const yaml = require("js-yaml");
const Path = require("path");
const mkdirp = require("mkdirp");
const stat = Promise.denodeify(require("fs").stat);
const writeFile = Promise.denodeify(require("fs").writeFile);
const waitOn = require("../shared/functions/utils").waitOn;

const absolutePath = (path) => {
  if (path.match(/^\.\./)) path = Path.join(__dirname, path);
  return Promise.resolve(path);
};

const fileExists = async (filename) => {
  let abs_path = await absolutePath(filename);
  return stat(abs_path)
    .then(() => {
      return true;
    })
    .catch(() => {
      return false;
    });
};

const writeJSON = async (filename, data) => {
  let abs_path = await absolutePath(filename);
  let match = abs_path.match(/^(.+)\/(.+)$/, "");
  let ready = await mkdirp(match[1]);
  let path = await waitOn(abs_path, ready);
  return writeFile(path, JSON.stringify(data))
    .then(() => {
      return true;
    })
    .catch((err) => {
      return err;
    });
};

const checkCompression = async (path) => {
  let abs_path = await absolutePath(path);
  let exists;
  if (!abs_path.endsWith(".gz")) {
    exists = await fileExists(abs_path);
  }
  if (exists) {
    return abs_path;
  } else {
    if (!abs_path.endsWith(".gz")) {
      abs_path += ".gz";
    }
    exists = await fileExists(abs_path);
    if (exists) {
      return abs_path;
    }
  }
  return undefined;
};

const readYaml = async (path) => {
  let abs_path = await absolutePath(path);
  let exists;
  if (!abs_path.endsWith(".gz")) {
    exists = await fileExists(abs_path);
  }
  if (exists) {
    let raw = fs.readFileSync(abs_path);
    return yaml.load(raw);
    // return read.sync(abs_path);
  } else {
    if (!abs_path.endsWith(".gz")) {
      abs_path += ".gz";
    }
    exists = await fileExists(abs_path);
    if (exists) {
      let raw = fs.readFileSync(abs_path);
      return yaml.load(zlib.gunzipSync(raw));
    }
  }
  return undefined;
};

const readYamlSync = (path) => {
  let raw = fs.readFileSync(path);
  if (path.endsWith(".gz")) {
    raw = zlib.gunzipSync(raw);
  }
  return yaml.load(raw);
};

const readJSON = readYaml;

// const readJSON = async (filename) => {
//   let abs_path = await absolutePath(filename);
//   return readFile(abs_path, "utf8")
//     .then(JSON.parse)
//     .catch((err) => {
//       return undefined;
//     });
// };

module.exports = {
  absolutePath,
  fileExists,
  readJSON,
  writeJSON,
  readYaml,
  readYamlSync,
  checkCompression,
};
