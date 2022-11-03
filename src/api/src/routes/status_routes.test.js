const app = require("../app");
const supertest = require("supertest");

test("GET /api/v1/status", async () => {
  const response = await supertest(app).get("/api/v1/status");
  expect(response.statusCode).toBe(200);
  expect(response.header["content-type"]).toBe(
    "application/json; charset=utf-8"
  );
  expect(response.body.status).toBe("LOADING");
});
