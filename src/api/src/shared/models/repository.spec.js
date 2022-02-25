const should = require('should');
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
});
