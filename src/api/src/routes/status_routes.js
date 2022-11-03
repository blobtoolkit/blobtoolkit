module.exports = function (app, db) {
  /**
   * @swagger
   * /api/v1/status:
   *   get:
   *     tags:
   *       - Status
   *     description: Returns API status and dataset count
   *     produces:
   *       - application/json
   *     responses:
   *       200:
   *         description: Check the API is running
   */
  app.get("/api/v1/status", async (req, res) => {
    res.setHeader("content-type", "application/json");
    const searchRoutes = require("./search_routes");
    let count = searchRoutes.total();
    let status = searchRoutes.status();
    if (count && typeof count === "number") {
      res.send({
        status,
        count,
      });
    } else {
      res.send({
        status,
      });
    }
  });
};
