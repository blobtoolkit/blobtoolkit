const appRoot = require('app-root-path');
const should = require('should');
const td = require('testdouble');
const config = require('../../config/main')
const Dataset = require('./dataset');

let dataset = undefined;
const id = 'ds1';
const blobDBFile = config.filePath + '/blobDB.json';
const filePath = config.filePath;
const blobDB = require(blobDBFile);

describe("Dataset model:", () => {
  beforeEach((done) => {
    dataset = new Dataset(id);
    done();
  })
  describe("loadBlobDB()", () => {
    it("should load a blobDB for a dataset", async () => {
      try {
        let blobDB = await dataset.loadBlobDB(blobDBFile);
        should.exist(blobDB);
      }
      catch(err) {
        should.not.exist(err);
      }
    });
  });
  describe("prepareMeta()", () => {
    beforeEach(()=>{
      dataset.blobDB = blobDB;
    })
    it("should convert a blobDB into a meta object", async () => {
      try {
        let meta = await dataset.prepareMeta(filePath);
        should.exist(meta);
        meta.should.have.property('filePath',filePath);
        meta.should.have.properties('id','name','records','record_type','fields');
        meta.id.should.be.String();
        meta.id.should.not.be.empty();
        meta.id.should.match(/^[\w\d]+$/);
        meta.name.should.be.String();
        meta.name.should.not.be.empty();
        meta.records.should.be.Number();
        meta.records.should.be.greaterThan(0);
        meta.record_type.should.be.String();
        meta.record_type.should.not.be.empty();
        meta.fields.should.be.Array();
        meta.fields.should.not.be.empty();
      }
      catch(err) {
        should.not.exist(err);
      }
    });
  });
  describe("storeMeta()", () => {
    let promise;
    beforeEach(()=>{
      dataset.blobDBFile = blobDBFile;
      promise = dataset.prepareMeta(filePath)
    })
    it("should store meta for a dataset", () => {
      promise.then(async () => {
        try {
          let success = await dataset.storeMeta();
          should.exist(success);
          success.should.be.true();
        }
        catch(err) {
          should.not.exist(err);
        }
      })
    });
  });
  describe("loadMeta()", () => {
    it("should load meta for a dataset", async () => {
      try {
        let meta = await dataset.loadMeta(id);
        should.exist(meta);
      }
      catch(err) {
        should.not.exist(err);
      }
    });
  });
  describe("storeLineages()", () => {
    let promise;
    beforeEach(() => {
      promise = dataset.loadMeta(id);
      dataset.blobDB = blobDB;
    });
    it("should store lineages for a dataset", () => {
      promise.then(async () => {
        try {
          let success = await dataset.storeLineages();
          should.exist(success);
          success.should.be.true();
        }
        catch(err) {
          console.log(err);
        }
      }).catch((err) => {
        console.log(err)
      });
    });
  });
  describe("addField()", () => {
    it("should add a field", () => {
      dataset.addField('a',{key:'value'});
      should.exist(dataset.fields);
      dataset.fields.should.be.a.Object();
      dataset.fields.should.have.property('a');
      dataset.fields['a'].should.be.a.Object();
      dataset.fields['a'].should.have.property('_id','a');
      dataset.fields['a'].should.have.property('_key','value');
    });
  });
  describe("addFields()", () => {
    it("should store fields for a dataset", () => {
      dataset.addFields([{id:'a',key:'value',children:[{id:'b'}]}]);
      should.exist(dataset.fields);
      dataset.fields.should.be.a.Object();
      dataset.fields.should.have.properties('b');
      dataset.fields['b'].should.be.a.Object();
      dataset.fields['b'].should.have.property('_id','b');
      dataset.fields['b'].should.have.property('_key','value');
    });
  });

});
