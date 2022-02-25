const express = require("express");
const path = require("path");
const PORT = process.env.BTK_PORT || process.env.BTK_CLIENT_PORT || "8080";
const BTK_API_PORT = process.env.BTK_API_PORT || "8000";
const BTK_API_URL =
  process.env.BTK_API_URL || `http://localhost:${BTK_API_PORT}/api/v1`;
const BTK_GA_ID = process.env.BTK_GA_ID || "UA-000000-01";
const BTK_GDPR_URL = process.env.BTK_GDPR_URL || undefined;
const app = express();

const ENV = {
  BTK_API_URL,
  BTK_GA_ID,
  BTK_GDPR_URL,
};

// set the view engine to ejs
app.set("views", __dirname + "/views");
app.set("view engine", "ejs");

// serve static assets normally
app.use("/view", express.static(path.resolve(__dirname, "public")));
app.use(
  "/manifest.json",
  express.static(path.resolve(__dirname, "public/manifest.json"))
);

// handle every other route with index.html, which will contain
// a script tag to your application's JavaScript file(s).
app.get("*", function (req, res) {
  // response.sendFile(path.resolve(__dirname, "view/index.html"));
  res.render("index", {
    variables: `window.process = ${JSON.stringify({ ENV })}`,
  });
});

app.listen(PORT);
console.log("server started on port " + PORT);
