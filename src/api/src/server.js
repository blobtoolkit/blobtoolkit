const OpenApiValidator = require("express-openapi-validator");
const swaggerUi = require("swagger-ui-express");
const YAML = require("yamljs");

const cookieParser = require("cookie-parser");

const express = require("express");
const compression = require("compression");
const bodyParser = require("body-parser");
const path = require("path");
const config = require("./config/main");
const logAccess = require("./functions/logger").logAccess;
const logError = require("./functions/logger").logError;
const resolve = require("path").resolve;

const apiSpec = path.join(__dirname, "api-spec.yaml");

let swaggerDocument = YAML.load(apiSpec);

swaggerDocument.servers[0].url = config.apiUrl;

const swaggerOptions = {
  customCss: ".swagger-ui .topbar { display: none }",
  customSiteTitle: `${config.title} API`,
};

const app = express();
app.use(compression());
if (config.cors) {
  const cors = require("cors");
  app.use(cors(config.cors));
}

app.use(express.urlencoded({ extended: true }));
app.use(express.text());
app.use(express.json());
app.use(cookieParser());

app.use((req, res, next) => {
  logAccess({ req });
  next();
});

// app.use(express.static(path.join(__dirname, "public")));
app.get("/api-spec", function (req, res) {
  res.header("Content-Type", "text/yaml");
  res.send(YAML.stringify(swaggerDocument, 8, 4));
});
if (process.pkg) {
  app.use(express.static(path.join(__dirname, "api-docs")));
  // serve swagger-ui.css directly as static path above only works for js
  app.get("/api-docs/swagger-ui.css", function (req, res) {
    res.header("Content-Type", "text/css");
    res.sendFile(path.join(__dirname, "api-docs", "swagger-ui.css"));
  });
}

app.locals.directory = resolve(config.filePath);

app.use(bodyParser.urlencoded({ extended: true }));

app.use(bodyParser.json());

const swaggerSetup = swaggerUi.setup(swaggerDocument, swaggerOptions);
app.get("/api-docs/index.html", swaggerSetup);
app.use("/api-docs", swaggerUi.serve);
app.get("/api-docs", swaggerSetup);

app.use(
  OpenApiValidator.middleware({
    apiSpec: swaggerDocument,
    validateRequests: {
      allowUnknownQueryParameters: true,
      removeAdditional: "failing",
    },
    validateResponses: false,
    operationHandlers: path.join(__dirname),
  })
);

app.use((err, req, res, next) => {
  let error = {
    message: err.message,
    errors: err.errors,
  };
  res.status(err.status || 500).json(error);
  logError({ ...error, req });
});

// require("./routes")(app, {});

module.exports = app;
