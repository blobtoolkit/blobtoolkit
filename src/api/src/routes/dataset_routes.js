const express = require("express");
const Repository = require("../models/repository");
const Dataset = require("../models/dataset");
const utils = require("../shared/functions/utils");

module.exports = function (app, db) {
  const datasetRoutes = express.Router();
  /**
   * @swagger
   * definitions:
   *   Dataset:
   *     properties:
   *       id:
   *         type: string
   *       name:
   *         type: string
   *       description:
   *         type: string
   *       records:
   *         type: integer
   *       record_type:
   *         type: string
   *       fields:
   *         type: array
   *         items:
   *           $ref: '#/definitions/Field'
   */
  /**
   * @swagger
   * parameters:
   *   dataset_id:
   *     in: path
   *     name: dataset_id
   *     type: string
   *     required: true
   *     description: Unique identifier of the requested dataset (e.g. ACVV01)
   */
  /**
   * @swagger
   * parameters:
   *   key:
   *     in: path
   *     name: key
   *     type: string
   *     required: true
   *     description: Dataset object key (e.g. assembly)
   */
  /**
   * @swagger
   * parameters:
   *   subkey:
   *     in: path
   *     name: subkey
   *     type: string
   *     required: true
   *     description: Key object subkey (e.g. accession)
   */
  /**
   * @swagger
   * /api/v1/dataset/all:
   *   get:
   *     tags:
   *       - Datasets
   *     description: Returns all datasets
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: An array of datasets
   *         schema:
   *           $ref: '#/definitions/Dataset'
   */
  datasetRoutes.get("/dataset/all", async (req, res) => {
    res.setHeader("content-type", "application/json");
    let repository = new Repository("default");
    let meta = await repository.loadMeta();
    res.json(meta);
  });
  /**
   * @swagger
   * /api/v1/dataset/id/{dataset_id}:
   *   get:
   *     tags:
   *       - Datasets
   *     description: Returns the dataset specified by dataset_id
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: A dataset object
   *         schema:
   *           $ref: '#/definitions/Dataset'
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   */
  datasetRoutes.get("/dataset/id/:dataset_id", async (req, res) => {
    res.setHeader("content-type", "application/json");
    let dataset = new Dataset(req.params.dataset_id);
    let meta = await dataset.loadMeta();
    if (meta.hasOwnProperty("id")) {
      res.json(meta);
    } else {
      res.sendStatus(404);
    }
  });
  /**
   * @swagger
   * /api/v1/dataset/id/{dataset_id}/{key}:
   *   get:
   *     tags:
   *       - Datasets
   *     description: Returns the dataset parameter specified by dataset_id and key
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: A dataset parameter
   *         schema:
   *           type: object
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   *       - $ref: "#/parameters/key"
   */
  datasetRoutes.get("/dataset/id/:dataset_id/:key", async (req, res) => {
    res.setHeader("content-type", "application/json");
    let dataset = new Dataset(req.params.dataset_id);
    let meta = await dataset.loadMeta();
    if (meta.hasOwnProperty(req.params.key)) {
      res.json(meta[req.params.key]);
    } else {
      let result = utils.nestedEntryByKeyValue(
        meta.fields,
        "_id",
        req.params.key,
        ["_data", "_children"]
      )[0];
      if (result) {
        res.json(result);
      } else {
        res.sendStatus(404);
      }
    }
  });

  /**
   * @swagger
   * /api/v1/dataset/id/{dataset_id}/{key}/{subkey}:
   *   get:
   *     tags:
   *       - Datasets
   *     description: Returns the dataset parameter specified by dataset_id, key and subkey
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: A dataset parameter
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   *       - $ref: "#/parameters/key"
   *       - $ref: "#/parameters/subkey"
   */
  datasetRoutes.get(
    "/dataset/id/:dataset_id/:key/:subkey",
    async (req, res) => {
      res.setHeader("content-type", "application/json");
      let dataset = new Dataset(req.params.dataset_id);
      let meta = await dataset.loadMeta();
      let json;
      if (meta.hasOwnProperty(req.params.key)) {
        json = meta[req.params.key];
      } else {
        let result = utils.nestedEntryByKeyValue(
          meta.fields,
          "_id",
          req.params.key,
          ["_data", "_children"]
        )[0];
        if (result) {
          json = result;
        } else {
          res.sendStatus(404);
        }
      }
      if (json.hasOwnProperty(req.params.subkey)) {
        res.json(json[req.params.subkey]);
      } else {
        res.sendStatus(404);
      }
    }
  );

  app.use("/api/v1", datasetRoutes);
};
