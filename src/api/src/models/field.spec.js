const should = require('should');
const td = require('testdouble');
const config = require('../../config/main')
const Field = require('./field');

describe("Field model:", () => {
  describe("loadData(id,index)", () => {
    it("should load all values for a field", async () => {
      let field = new Field('gc',{id:'ds1'});
      let result = await field.loadData();
      result.should.be.a.Object();
      result.should.have.property('values');
    });
  });
  describe("loadDataAtIndex()", () => {
    it("should load a specific value for a field", async () => {
      let field = new Field('gc',{id:'ds1'});
      let result = await field.loadDataAtIndex(7);
      result.should.be.a.Object();
      result.should.have.property('values');
      result.values.length.should.equal(1);
      result.values[0].should.equal(0.2801);
    });
    it("should load a range of values for a field", async () => {
      let field = new Field('gc',{id:'ds1'});
      let result = await field.loadDataAtIndex('5-7');
      result.should.be.a.Object();
      result.should.have.property('values');
      result.values.length.should.equal(3);
      result.values[0].should.equal(0.1944);
    });
    it("should load comma-separated ranges", async () => {
      let field = new Field('gc',{id:'ds1'});
      let result = await field.loadDataAtIndex('1,5-7');
      result.should.be.a.Object();
      result.should.have.property('values');
      result.values.length.should.equal(4);
      result.values[0].should.equal(0.2623);
      result.values[3].should.equal(0.2801);
    });
    it("should reject incomplete ranges", async () => {
      let field = new Field('gc',{id:'ds1'});
      let result = await field.loadDataAtIndex('1,-7');
      should.not.exist(result);
    });
    it("should reject reversed ranges", async () => {
      let field = new Field('gc',{id:'ds1'});
      let result = await field.loadDataAtIndex('8-7');
      should.not.exist(result);
    });
    it("should return translated value for key-value fields", async () => {
      let field = new Field('bestsum_family',{id:'ds1'});
      let result = await field.loadDataAtIndex('5');
      result.should.be.a.Object();
      result.should.have.property('values');
      result.values.length.should.equal(1);
      result.values[0].should.match('Hypsibiidae');
    });
  });
  describe("scaleType(value)", (done) => {
    let field;
    beforeEach((done) => {
      field = new Field('gc',{id:'ds1'},{range:[1,4]});
      done();
    })
    it("should return the name of the scale function if called with no value", (done) => {
      let result = field.scaleType();
      result.should.be.a.String();
      result.should.match('scaleLinear');
      done();
    });
    it("should change the scale function if a valid value is passed", (done) => {
      let result = field.scaleType('scaleLog');
      result.should.be.a.String();
      result.should.match('scaleLog');
      result = field.scaleType('scaleSqrt');
      result.should.be.a.String();
      result.should.match('scaleSqrt');
      done();
    });
    it("should return undefined if an invalid value is passed", (done) => {
      let result = field.scaleType('scaleInvalid');
      should.not.exist(result);
      done();
    });
  });
  describe("range(arr)", () => {
    let field;
    beforeEach((done) => {
      field = new Field('gc',{id:'ds1'},{range:[0,2]});
      done();
    })
    it("should set range limits for the field", (done) => {
      field._range.should.be.a.Array();
      field._range.should.match([0,2]);
      done();
    });
    it("should return range when called with no values", (done) => {
      let range = field.range();
      range.should.be.a.Array();
      range.should.match([0,2]);
      done();
    });
    it("should return range when called with an empty array", (done) => {
      let range = field.range([]);
      range.should.be.a.Array();
      range.should.match([0,2]);
      done();
    });
    it("should set allow inverted max and min", (done) => {
      field.range([3,1]);
      field._range.should.be.a.Array();
      field._range.should.match([1,3]);
      done();
    });
    it("should find max and min if array is too long", (done) => {
      field.range([2,3,0,1]);
      field._range.should.be.a.Array();
      field._range.should.match([0,3]);
      done();
    });
    it("should set max and min to a single value if only one value is passed", (done) => {
      field.range([2]);
      field._range.should.be.a.Array();
      field._range.should.match([2,2]);
      done();
    });
  });
  describe("rangeHigh(value)", (done) => {
    let field;
    beforeEach((done) => {
      field = new Field('gc',{id:'ds1'},{range:[1,3]});
      done();
    })
    it("should increase upper range limit for the field", (done) => {
      field.rangeHigh(4);
      field._range.should.be.a.Array();
      field._range.should.match([1,4]);
      done();
    });
    it("should reduce upper range limit for the field", (done) => {
      field.rangeHigh(2);
      field._range.should.be.a.Array();
      field._range.should.match([1,2]);
      done();
    });
    it("should reduce the lower range limit if below existing range", (done) => {
      field.rangeHigh(0);
      field._range.should.be.a.Array();
      field._range.should.match([0,0]);
      done();
    });
    it("should return upper range limit when called with no values", (done) => {
      let range = field.rangeHigh();
      range.should.be.a.Number();
      range.should.equal(3);
      done();
    });
    it("should return upper range limit when called with an invalid datatype", (done) => {
      let range = field.rangeHigh('2');
      range.should.be.a.Number();
      range.should.equal(3);
      done();
    });
  });
  describe("rangeLow(value)", (done) => {
    let field;
    beforeEach((done) => {
      field = new Field('gc',{id:'ds1'},{range:[1,3]});
      done();
    })
    it("should reduce lower range limit for the field", (done) => {
      field.rangeLow(0);
      field._range.should.be.a.Array();
      field._range.should.match([0,3]);
      done();
    });
    it("should increase lower range limit for the field", (done) => {
      field.rangeLow(2);
      field._range.should.be.a.Array();
      field._range.should.match([2,3]);
      done();
    });
    it("should increase the upper range limit if value is above existing range", (done) => {
      field.rangeLow(4);
      field._range.should.be.a.Array();
      field._range.should.match([4,4]);
      done();
    });
    it("should return lower range limit when called with no values", (done) => {
      let range = field.rangeLow();
      range.should.be.a.Number();
      range.should.equal(1);
      done();
    });
    it("should return lower range limit when called with an invalid datatype", (done) => {
      let range = field.rangeLow('2');
      range.should.be.a.Number();
      range.should.equal(1);
      done();
    });
  });
  describe("FilterToList(arr)", () => {
    let field;
    beforeEach(() => {
      field = new Field('gc',{id:'ds1'});
      field.filters['default']._range = [0.3,0.7];
    });
    it("should return an array of indices for values that pass filter", async () => {
      let result = await field.filterToList();
      result.should.be.a.Array();
      result.length.should.equal(3);
      result.should.match([3,4,9]);
    });
    it("should return inverse if inclusive is false", async () => {
      field.filters['default'].inclusive(false);
      let result = await field.filterToList();
      result.should.be.a.Array();
      result.length.should.equal(7);
      result.should.match([0,1,2,5,6,7,8]);
    });
    it("should restrict results to passed list", async () => {
      let result = await field.filterToList([1,2,3,4]);
      result.should.be.a.Array();
      result.length.should.equal(2);
      result.should.match([3,4]);
    });
  });
  describe("applyListFilter(arr)", () => {
    it("should return an array values", async () => {
      let field = new Field('gc',{id:'ds1'});
      let list = [1,3,5,7];
      let result = await field.applyListFilter(list);
      result.should.be.a.Array();
      result.length.should.equal(4);
      result.should.match([0.2623,0.3155,0.1944,0.2801]);
    });
    it("should return values from key value pairs", async () => {
      let field = new Field('bestsum_family',{id:'ds1'});
      let list = [1,3,5,7];
      let result = await field.applyListFilter(list);
      result.should.be.a.Array();
      result.length.should.equal(4);
      result.should.match(['Microbacteriaceae','Microbacteriaceae','Hypsibiidae','unresolved']);
    });
  });

});
