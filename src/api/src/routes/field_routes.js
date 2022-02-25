const appRoot = require("app-root-path");
const Field = require("../models/field");
const utils = require("../shared/functions/utils");

module.exports = function (app, db) {
  /**
   * @swagger
   * definitions:
   *   Field:
   *     properties:
   *       id:
   *         type: string
   *         description: a unique identifier
   *         required: true
   *       name:
   *         description: a human-readable name
   *         type: string
   *         required: true
   *       description:
   *         description: a brief description
   *         type: string
   *         required: false
   *       datatype:
   *         description: the datatype must be one of integer, float, boolean or string
   *         type: string
   *         enum:
   *           - float
   *           - integer
   *           - mixed
   *           - string
   *         required: true
   *       type:
   *         description: the type must be one of variable, category, label or identifier
   *         type: string
   *         enum:
   *           - identifier
   *           - variable
   *           - category
   *           - array
   *           - multiarray
   *         required: true
   *       range:
   *         description: (optional) an array of length 2 containing minimum and maximum values
   *         type: array
   *         items:
   *           type: number
   *           format: float
   *         minItems: 2
   *         maxItems: 2
   *         required: false
   *       preload:
   *         description: (optional) a flag to indicate whether this field should be loaded for visualisation without user prompting
   *         type: boolean
   *         required: false
   *       active:
   *         description: (optional) a flag to indicate whether this field should be part of preview set
   *         type: boolean
   *         required: false
   *       scale:
   *         description: (optional) scale to use when binning/presenting data
   *         type: string
   *         enum:
   *           - scaleLinear
   *           - scaleLog
   *           - scaleSqrt
   *         required: false
   *       clamp:
   *         description: (optional) clamp low values to number (used to allow use of log scale when range contains zero)
   *         type: [number, boolean]
   *         required: false
   *       linked_field:
   *         description: (optional) field ID to which the current field is linked
   *         type: string
   *         pattern: ^[A-Za-z0-9_.-]+$
   *         required: false
   *       children:
   *         description: (optional) an array of nested fields
   *         type: array
   *         items:
   *           type: string
   *         required: false
   *       data:
   *         description: (optional) an array of nested data
   *         type: array
   *         items:
   *           type: object
   *         required: false
   */
  /**
   * @swagger
   * definitions:
   *   Value:
   *     properties:
   *       values:
   *         type: array
   *         description: an array of values
   *         required: true
   *       keys:
   *         type: array
   *         items:
   *           type: object
   *         description: an array of key-value pairs
   *         required: false
   */
  /**
   * @swagger
   * parameters:
   *   field_id:
   *     in: path
   *     name: field_id
   *     type: string
   *     required: true
   *     description: Unique identifier of the requested field (e.g. gc)
   */
  /**
   * @swagger
   * parameters:
   *   index:
   *     in: path
   *     name: index
   *     type: string
   *     pattern: ^\d+[,-\d]*
   *     required: true
   *     description: index of the requested record(s) (e.g. 72 or 11,56,99-101)
   */
  /**
   * @swagger
   * /api/v1/field/{dataset_id}/{field_id}:
   *   get:
   *     tags:
   *       - Fields
   *     description: Returns all values for a field
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: An object containing an array of variable values
   *         schema:
   *           $ref: '#/definitions/Value'
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   *       - $ref: "#/parameters/field_id"
   */
  app.get("/api/v1/field/:dataset_id/:field_id", async (req, res) => {
    res.setHeader("content-type", "application/json");
    //res.setHeader('Access-Control-Allow-Origin', 'http://localhost:8080');
    let field = new Field(req.params.field_id, { id: req.params.dataset_id });
    let data = await field.loadData();
    if (data && data.hasOwnProperty("values")) {
      res.json(data);
    } else {
      res.sendStatus(404);
    }
  });
  /**
   * @swagger
   * /api/v1/field/{dataset_id}/{field_id}/{index}:
   *   get:
   *     tags:
   *       - Fields
   *     description: Returns specific values for a field
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: An object containing an array of variable values
   *         schema:
   *           $ref: '#/definitions/Value'
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   *       - $ref: "#/parameters/field_id"
   *       - $ref: "#/parameters/index"
   */
  app.get("/api/v1/field/:dataset_id/:field_id/:index", async (req, res) => {
    res.setHeader("content-type", "application/json");
    let field = new Field(req.params.field_id, { id: req.params.dataset_id });
    let data = await field.loadDataAtIndex(req.params.index);
    if (data) {
      res.json(data);
    } else {
      res.sendStatus(404);
    }
  });
};
