const express = require("express");
const compression = require("compression");
const bodyParser = require("body-parser");
const path = require("path");
const config = require("./config/main");
const logAccess = require("./functions/logger").logAccess;
const logError = require("./functions/logger").logError;
const resolve = require("path").resolve;
const app = express();
const searchRoutes = require("./routes/search_routes");

app.use(compression());
app.use(express.static(path.join(__dirname, "public")));
if (process.pkg) {
  app.use(express.static(path.join(__dirname, "api-docs")));
}
app.locals.directory = resolve(config.filePath);
if (config.cors) {
  const cors = require("express-cors");
  app.use(cors(config.cors));
}
app.use(bodyParser.urlencoded({ extended: true }));

app.use(bodyParser.json());

app.use((req, res, next) => {
  logAccess({ req });
  next();
});

app.use((err, req, res, next) => {
  let error = {
    message: err.message,
    errors: err.errors,
  };
  res.status(err.status || 500).json(error);
  logError({ ...error, req });
});

require("./routes")(app, {});

module.exports = app;
