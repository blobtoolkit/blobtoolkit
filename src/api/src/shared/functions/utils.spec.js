const should = require('should');
const utils = require('./utils');
const Promise = require('promise');
const td = require('testdouble');

describe("utils:", () => {
  describe("waitOn(value,promise):", () => {
    it("should wait on a promise", async () => {
      let delay = (time) => {
        return new Promise((fulfill) => {
          setTimeout(fulfill, time);
        });
      }
      let promise = await delay(10);
      return utils.waitOn('done',promise)
        .then((value) => {
          should.exist(value);
          value.should.match('done')
        });
    });
  });
  describe("asKeyValue(arr):", () => {
    it("should convert an array with repeated elements to an object with keys and values", async () => {
      let obj = await utils.asKeyValue(['a','b','a','c'])
      should.exist(obj);
      obj.should.have.properties('keys','values');
      obj.should.have.property('keys',['a','b','c']);
      obj.should.have.property('values',[0,1,0,2]);
    });
  });
  describe("groupValuesBy(arr,by,val):", () => {
    it("should group an array of objects 'by' object key, returning arrays of values for 'val'", async () => {
      let obj = await utils.groupValuesBy([{x:'A',y:3},{x:'B',y:2},{x:'A',y:1}],'x','y')
      should.exist(obj);
      obj.should.have.properties('A','B');
      obj.A.should.have.length(2);
      obj.A.should.match([3,1])
      obj.B.should.have.length(1);
      obj.B.should.match([2])
    });
  });
  describe("valueAtPath(obj,path):", () => {
    it("should return the value from a nested object following an array of keys in 'path'", async () => {
      let val = await utils.valueAtPath({id:'A',data:{coords:{x:3,y:4}}},['data','coords','y'])
      should.exist(val);
      val.should.be.a.Number();
      val.should.equal(4)
    });
  });
  describe("newField(id,options):", () => {
    it("should create an object with id and optional properties", async () => {
      let field = await utils.newField('x',{description:'test'})
      should.exist(field);
      field.should.have.property('id','x')
      field.should.have.property('name','x')
      field.should.have.property('description','test')
    });
  });
  describe("entryByKeyValue(arr,key,value):", () => {
    let arr = [{id:1,type:'a'},{id:2,type:'b'},{id:3,type:'a'}];
    it("should return a single object with unique key value", (done) => {
      let result = utils.entryByKeyValue(arr,'type','b')
      should.exist(result);
      result.should.be.a.Array();
      result.should.have.length(1);
      result[0].should.have.property('id',2)
      done();
    });
    it("should return an array of matching objects", (done) => {
      let result = utils.entryByKeyValue(arr,'type','a')
      should.exist(result);
      result.should.be.a.Array();
      result.should.have.length(2);
      result[0].should.be.a.Object();
      done();
    });

    it("should return undefined if no matching object is found", (done) => {
      let result = utils.entryByKeyValue(arr,'type','c')
      should.not.exist(result);
      done();
    });
  });
  describe("nestedEntryByKeyValue(arr,key,value,nestarr):", () => {
    let arr = [{id:1,type:'a'},{id:2,type:'b',children:[{id:4,type:'d'},{id:5,type:'d',data:[{id:8,type:'e'}]}]},{id:3,type:'a',data:[{id:6,type:'e'},{id:7,type:'d'}]}];
    let nestarr = ['data','children'];
    let entryByKeyValue;
    it("should return a top-level object", (done) => {
      let result = utils.nestedEntryByKeyValue(arr,'type','b',nestarr);
      should.exist(result);
      result.should.be.a.Array();
      result.should.have.length(1);
      result[0].should.have.property('id',2);
      done();
    });
    it("should return a nested object", (done) => {
      let result = utils.nestedEntryByKeyValue(arr,'id','4',nestarr);
      should.exist(result);
      result.should.be.a.Array();
      result.should.have.length(1);
      result[0].should.have.property('type','d');
      done();
    });
    it("should return a deeply nested object", (done) => {
      let result = utils.nestedEntryByKeyValue(arr,'id','8',nestarr);
      should.exist(result);
      result.should.be.a.Array();
      result.should.have.length(1);
      result[0].should.have.property('type','e');
      done();
    });
    it("should return a objects from different nesting levels", (done) => {
      let result = utils.nestedEntryByKeyValue(arr,'type','e',nestarr);
      should.exist(result);
      result.should.be.a.Array();
      result.should.have.length(2);
      result[0].should.have.property('id',8);
      result[1].should.have.property('id',6);
      done();
    });
  });
});
