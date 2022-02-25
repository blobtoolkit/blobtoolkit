const express = require("express");
const path = require("path");
const swaggerJSDoc = require("swagger-jsdoc");
const swaggerUI = require("swagger-ui-express");
const config_swagger = require("../config/swagger");
const YAML = require("yamljs");

const swaggerSpec = swaggerJSDoc(config_swagger.options);

module.exports = function (app, db) {
  // serve swagger

  if (process.pkg) {
    app.get("/api-docs/swagger-ui.css", function (req, res) {
      res.header("Content-Type", "text/css");
      res.sendFile(
        path.resolve(path.join(__dirname, "..", "api-docs", "swagger-ui.css"))
      );
    });
  }
  app.use("/api-docs", swaggerUI.serve, swaggerUI.setup(swaggerSpec));
  app.get("/swagger.json", function (req, res) {
    res.setHeader("Content-Type", "application/json");
    res.send(swaggerSpec);
  });
  app.get("/api-spec", function (req, res) {
    res.header("Content-Type", "text/yaml");
    res.send(YAML.stringify(swaggerSpec, 8, 4));
  });
};
