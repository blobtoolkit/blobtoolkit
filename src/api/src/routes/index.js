const apiRoutes = require("./api_routes");
const datasetRoutes = require("./dataset_routes");
const fieldRoutes = require("./field_routes");
const identifierRoutes = require("./identifier_routes");
const statusRoutes = require("./status_routes");
const imageRoutes = require("./image_routes");
const searchRoutes = require("./search_routes");
// const sliceRoutes = require('./slice_routes');
const summaryRoutes = require("./summary_routes");
const swaggerRoutes = require("./swagger_routes");

module.exports = async function (app, db) {
  apiRoutes(app, db);
  datasetRoutes(app, db);
  fieldRoutes(app, db);
  statusRoutes(app, db);
  identifierRoutes(app, db);
  imageRoutes(app, db);
  searchRoutes.routes(app, db);
  // sliceRoutes(app, db);
  summaryRoutes(app, db);
  swaggerRoutes(app, db);
};
