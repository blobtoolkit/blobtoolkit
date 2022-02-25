const express = require('express');

module.exports = function(app) {
  const apiRoutes = express.Router();
  // Set url for API group routes
  app.use('/api/v1', apiRoutes);
};
