const should = require('should');
const main_config = require('../../config/main')
const Repository = require('./repository');

let repository = function(){};
const id = 'default';

describe("Repository:", () => {
  beforeEach((done) => {
    repository = new Repository(id);
    done();
  })
  describe("getId()", () => {
    it("should show the repository ID", () => {
      value = repository.getId();
      value.should.match(id);
    });
  });
  describe("loadMeta()", () => {
    it("should load metadata for a repository", async () => {
      meta = await repository.loadMeta();
      meta[0].id.should.match('ds1');
    });
  });
});
