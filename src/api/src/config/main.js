require("dotenv").config();

const appRoot = require("app-root-path");
const BTK_HOST = process.env.BTK_HOST || "localhost";
const BTK_CLIENT_PORT =
  Number(process.env.BTK_CLIENT_PORT) || Number(process.env.BTK_PORT) || 8080;
const BTK_API_PORT = Number(process.env.BTK_API_PORT) || 8000;
const BTK_API_URL =
  process.env.BTK_API_URL || BTK_HOST + ":" + BTK_API_PORT + "/api/v1";
const BTK_HTTPS = String(process.env.BTK_HTTPS) === "true";
const BTK_ORIGINS = process.env.BTK_ORIGINS
  ? process.env.BTK_ORIGINS.split(" ")
  : [
      "localhost",
      "null",
      BTK_HOST,
      (BTK_HTTPS ? "https" : "http") + "://" + BTK_HOST + ":" + BTK_CLIENT_PORT,
    ];
const FILE_PATH = process.env.BTK_FILE_PATH || appRoot + "/demo";

module.exports = {
  // setting port for server
  client_port: BTK_CLIENT_PORT,
  // setting port for server
  api_port: BTK_API_PORT,
  // flag to use https
  https: BTK_HTTPS,
  // Cors settings
  cors: {
    allowedOrigins: BTK_ORIGINS,
  },
  // API URL
  apiUrl: BTK_API_URL,
  // url basename
  basename: process.env.BTK_BASENAME || "",
  // search index reload key
  reloadKey: process.env.BTK_RELOAD_KEY || "",
  // path to read flatfiles
  filePath: FILE_PATH,
  // Path to write flatfiles
  outFilePath: process.env.BTK_OUT_FILE_PATH || appRoot + "/files/out",
  // version
  version: process.env.BTK_VERSION || "v3.0.0",
  // hostname
  hostname: BTK_HOST,
  // API URL
  apiUrl:
    process.env.BTK_API_URL ||
    (BTK_HTTPS ? "https" : "http") +
      "://" +
      BTK_HOST +
      ":" +
      BTK_API_PORT +
      "/api/v1",
  mode: process.env.NODE_ENV || "test",
  // static threshold
  staticThreshold: process.env.BTK_STATIC_THRESHOLD || 100000,
  // nohit threshold
  nohitThreshold: process.env.BTK_NOHIT_THRESHOLD || 1000000,
  // circle limit
  circleLimit: process.env.BTK_CIRCLE_LIMIT || 100000,
  // SSL
  keyFile: process.env.BTK_KEYFILE || "",
  certFile: process.env.BTK_CERTFILE || "",
  ga_id: process.env.BTK_GA_ID || "UA-000000-01",
  gdpr_url: process.env.BTK_GDPR_URL,
  dataset_table: process.env.BTK_DATASET_TABLE || false,
  target_tree: process.env.BTK_TARGET_TREE || false,
  message: process.env.BTK_MESSAGE || false,
  use_default: process.env.BTK_USE_DEFAULT_LINKS || false,
  disableHostCheck: String(process.env.BTK_DISABLE_HOST_CHECK) === "true",
  aboutUrl: process.env.BTK_ABOUT_URL || "http://blobtoolkit.genomehubs.org",
};
