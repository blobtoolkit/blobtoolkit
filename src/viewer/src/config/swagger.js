const main = require('./main')
const appRoot = require('app-root-path');

// swagger definition
const swaggerDefinition = {
  info: {
    title: 'BlobToolKit API',
    version: main.version,
    description: 'A RESTful API for BlobToolKit',
  },
  host: main.hostname + (main.hostname == 'localhost' ? ':' + main.api_port : ''),
  basePath: '/',
};

// options for the swagger docs
const options = {
  // import swaggerDefinitions
  swaggerDefinition: swaggerDefinition,
  // path to the API docs
  apis: [appRoot + '/src/server/routes/*.js'],
};

module.exports = {
  options: options
}
