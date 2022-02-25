const checkCompression = require("../functions/io").checkCompression;

//const sharp = require("sharp");

// const resize = async (file, width, height) => {
//   height = height || width;
//   let buffer = await sharp(file)
//     .resize(width, height, {
//       fit: sharp.fit.inside,
//       withoutEnlargement: true,
//     })
//     .toFormat("png")
//     .toBuffer()
//     .then((outputBuffer) => outputBuffer);
//   return buffer;
// };

/**
 * @swagger
 * parameters:
 *   view:
 *     in: path
 *     name: view
 *     type: string
 *     enum:
 *       - blob
 *       - cumulative
 *       - kite
 *       - snail
 *     required: true
 *     description: view type
 */
/**
 * @swagger
 * parameters:
 *   type:
 *     in: path
 *     name: type
 *     type: string
 *     enum:
 *       - circle
 *       - hex
 *       - kite
 *       - none
 *       - square
 *     required: true
 *     description: blob view plot type
 */
/**
 * @swagger
 * parameters:
 *   format:
 *     in: query
 *     name: format
 *     type: string
 *     enum:
 *       - png
 *       - svg
 *     required: false
 *     description: image format
 *     default: png
 */
/**
 * @swagger
 * parameters:
 *   width:
 *     in: query
 *     name: width
 *     type: integer
 *     required: false
 *     description: png image width
 */

module.exports = function (app, db) {
  var directory = app.locals.directory;

  /**
   * @swagger
   * /api/v1/image/{dataset_id}:
   *   get:
   *     tags:
   *       - Images
   *     description: Returns an image of a dataset view
   *     produces:
   *       - image/png
   *       - image/svg+xml
   *     responses:
   *       200:
   *         description: OK
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   */
  app.get("/api/v1/image/:dataset_id", async (req, res) => {
    file = await checkCompression(
      directory + "/" + req.params.dataset_id + "/blob.circle.png"
    );
    if (file) {
      res.setHeader("content-type", "image/png");
      if (file.endsWith(".gz")) {
        res.setHeader("content-encoding", "gzip");
      }
      res.sendFile(file);
    } else {
      res.sendStatus(404);
    }
  });

  /**
   * @swagger
   * /api/v1/image/{dataset_id}/{view}:
   *   get:
   *     tags:
   *       - Images
   *     description: Returns an image of a dataset view
   *     produces:
   *       - image/png
   *       - image/svg+xml
   *     responses:
   *       200:
   *         description: OK
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   *       - $ref: "#/parameters/view"
   *       - $ref: "#/parameters/format"
   *       - $ref: "#/parameters/width"
   */
  app.get("/api/v1/image/:dataset_id/:view", async (req, res) => {
    let type = "";
    if (req.params.view == "blob") {
      type = ".circle";
    }
    if (req.query.format == "svg") {
      let file =
        directory +
        "/" +
        req.params.dataset_id +
        "/" +
        req.params.view +
        type +
        ".svg";
      file = await checkCompression(file);
      if (file) {
        res.setHeader("content-type", "image/svg+xml");
        if (file.endsWith(".gz")) {
          res.setHeader("content-encoding", "gzip");
        }
        res.sendFile(file);
      } else {
        res.sendStatus(404);
      }
    } else {
      res.setHeader("content-type", "image/png");
      let file =
        directory +
        "/" +
        req.params.dataset_id +
        "/" +
        req.params.view +
        type +
        ".png";
      file = await checkCompression(file);
      if (file) {
        res.setHeader("content-type", "image/png");
        if (file.endsWith(".gz")) {
          res.setHeader("content-encoding", "gzip");
        }
        res.sendFile(file);
      } else {
        res.sendStatus(404);
      }
    }
  });

  /**
   * @swagger
   * /api/v1/image/{dataset_id}/blob/{type}:
   *   get:
   *     tags:
   *       - Images
   *     description: Returns an image of a blob view plot type
   *     produces:
   *       - image/png
   *       - image/svg+xml
   *     responses:
   *       200:
   *         description: OK
   *     parameters:
   *       - $ref: "#/parameters/dataset_id"
   *       - $ref: "#/parameters/type"
   *       - $ref: "#/parameters/format"
   *       - $ref: "#/parameters/width"
   */
  app.get("/api/v1/image/:dataset_id/blob/:type", async (req, res) => {
    let type = req.params.type;
    if (type == "none") type = "circle";
    if (req.query.format == "svg") {
      let file =
        directory + "/" + req.params.dataset_id + "/blob." + type + ".svg";
      file = await checkCompression(file);
      if (file) {
        res.setHeader("content-type", "image/svg+xml");
        if (file.endsWith(".gz")) {
          res.setHeader("content-encoding", "gzip");
        }
        res.sendFile(file);
      } else {
        res.sendStatus(404);
      }
    } else {
      let file =
        directory + "/" + req.params.dataset_id + "/blob." + type + ".png";
      file = await checkCompression(file);
      if (file) {
        res.setHeader("content-type", "image/png");
        if (file.endsWith(".gz")) {
          res.setHeader("content-encoding", "gzip");
        }
        res.sendFile(file);
      } else {
        res.sendStatus(404);
      }
    }
  });
};
