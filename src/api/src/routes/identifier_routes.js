const checkCompression = require("../functions/io").checkCompression;

module.exports = function (app, db) {
  var directory = app.locals.directory;

  /**
   * @swagger
   * /api/v1/identifiers/{dataset_id}:
   *   get:
   *     tags:
   *       - Identifiers
   *     description: Returns all identifiers for a dataset
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: An object containing an array of variable values
   *         schema:
   *           $ref: '#/definitions/Value'
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   */
  app.get("/api/v1/identifiers/:dataset_id", async (req, res) => {
    file = await checkCompression(
      directory + "/" + req.params.dataset_id + "/identifiers.json"
    );
    if (file) {
      res.setHeader("content-type", "application/json");
      if (file.endsWith(".gz")) {
        res.setHeader("content-encoding", "gzip");
      }
      res.sendFile(file);
    } else {
      res.sendStatus(404);
    }
  });
};
