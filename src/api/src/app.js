const config = require("./config/main");
const app = require("./server");
const searchRoutes = require("./routes/search_routes");

if (config.https) {
  const https = require("https");
  const fs = require("fs");
  var options = {
    key: fs.readFileSync(config.keyFile),
    cert: fs.readFileSync(config.certFile),
  };
  var secureServer = https
    .createServer(options, app)
    .listen(config.api_port, async function () {
      console.log("blobtoolkit-api started on https port: " + config.api_port);
      await searchRoutes.loadIndex();
    });
} else {
  app.listen(config.api_port, async function () {
    console.log("blobtoolkit-api started on http port: " + config.api_port);
    await searchRoutes.loadIndex();
  });
}
module.exports = app;
