const getStatus = async (req, res) => {
  res.setHeader("content-type", "application/json");
  const search = require("./search");
  let count = search.total();
  let status = search.status();
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
};

module.exports = {
  getStatus,
};
