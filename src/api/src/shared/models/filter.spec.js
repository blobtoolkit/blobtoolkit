const should = require('should');
const config = require('../../config/main')
const Filter = require('./filter');
const Field = require('./field');

describe("Filter model:", () => {
/*  describe("scaleType(value)", (done) => {
    let filter;
    beforeEach((done) => {
      filter = new Filter('gc','ds1');
      done();
    })
    it("should return the name of the scale function if called with no value", (done) => {
      let result = filter.scaleType();
      result.should.be.a.String();
      result.should.match('scaleLinear');
      done();
    });
    it("should change the scale function if a valid value is passed", (done) => {
      let result = filter.scaleType('scaleLog');
      result.should.be.a.String();
      result.should.match('scaleLog');
      result = filter.scaleType('scaleSqrt');
      result.should.be.a.String();
      result.should.match('scaleSqrt');
      done();
    });
    it("should return undefined if an invalid value is passed", (done) => {
      let result = filter.scaleType('scaleInvalid');
      should.not.exist(result);
      done();
    });

  });*/
/*  describe("limits(arr)", () => {
    let filter;
    beforeEach((done) => {
      filter = new Filter('gc','ds1');
      done();
    })
    it("should set range limits for the filter", (done) => {
      filter.limits([0,2]);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([0,2]);
      done();
    });
    it("should return limits when called with no values", (done) => {
      filter.limits([0,2]);
      let limits = filter.limits();
      limits.should.be.a.Array();
      limits.should.match([0,2]);
      done();
    });
    it("should return limits when called with an empty array", (done) => {
      filter.limits([0,2]);
      let limits = filter.limits([]);
      limits.should.be.a.Array();
      limits.should.match([0,2]);
      done();
    });
    it("should set allow inverted max and min", (done) => {
      filter.limits([3,1]);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([1,3]);
      done();
    });
    it("should find max and min if array is too long", (done) => {
      filter.limits([2,3,0,1]);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([0,3]);
      done();
    });
    it("should set max and min to a single value if only one value is passed", (done) => {
      filter.limits([2]);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([2,2]);
      done();
    });
    it("should set range to match limits", (done) => {
      filter.limits([2,3,1]);
      filter._range.should.be.a.Array();
      filter._range.should.match([1,3]);
      done();
    });
  });
  describe("limitsHigh(value)", (done) => {
    let filter;
    beforeEach((done) => {
      filter = new Filter('gc','ds1');
      filter.limits([1,3]);
      done();
    })
    it("should increase upper limit for the filter", (done) => {
      filter.limitsHigh(4);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([1,4]);
      done();
    });
    it("should reduce upper limit for the filter", (done) => {
      filter.limitsHigh(2);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([1,2]);
      done();
    });
    it("should reduce the lower limit if below existing range", (done) => {
      filter.limitsHigh(0);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([0,0]);
      done();
    });
    it("should return upper limit when called with no values", (done) => {
      let limit = filter.limitsHigh();
      limit.should.be.a.Number();
      limit.should.equal(3);
      done();
    });
    it("should return upper limit when called with an invalid datatype", (done) => {
      let limit = filter.limitsHigh('2');
      limit.should.be.a.Number();
      limit.should.equal(3);
      done();
    });
    it("should update range if upper limit is changed", (done) => {
      filter.limitsHigh(2);
      filter._range.should.be.a.Array();
      filter._range.should.match([1,2]);
      done();
    });
  });
  describe("limitsLow(value)", (done) => {
    let filter;
    beforeEach((done) => {
      filter = new Filter('gc','ds1');
      filter.limits([1,3]);
      done();
    })
    it("should reduce lower limit for the filter", (done) => {
      filter.limitsLow(0);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([0,3]);
      done();
    });
    it("should increase lower limit for the filter", (done) => {
      filter.limitsLow(2);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([2,3]);
      done();
    });
    it("should increase the upper limit if value is above existing range", (done) => {
      filter.limitsLow(4);
      filter._limits.should.be.a.Array();
      filter._limits.should.match([4,4]);
      done();
    });
    it("should return lower limit when called with no values", (done) => {
      let limit = filter.limitsLow();
      limit.should.be.a.Number();
      limit.should.equal(1);
      done();
    });
    it("should return upper limit when called with an invalid datatype", (done) => {
      let limit = filter.limitsLow('2');
      limit.should.be.a.Number();
      limit.should.equal(1);
      done();
    });
    it("should update range if range becomes out of bounds", (done) => {
      filter.limitsLow(2);
      filter._range.should.be.a.Array();
      filter._range.should.match([2,3]);
      done();
    });
  });*/
  describe("range(arr)", (done) => {
    let field;
    let filter;
    beforeEach((done) => {
      field = new Field('gc','ds1');
      filter = field.filters['default'];
      done();
    })
    it("should set range limits for the filter", (done) => {
      field._range = [0,4];
      filter.range([1,3,2]);
      filter._range.should.be.a.Array();
      filter._range.should.match([1,3]);
      done();
    });
    it("should keep range within limits", (done) => {
      field._range = [1,2];
      filter.range([2,3,1,0]);
      filter._range.should.be.a.Array();
      filter._range.should.match([1,2]);
      done();
    });
  });
  describe("rangeHigh(value)", (done) => {
    let field;
    let filter;
    beforeEach((done) => {
      field = new Field('gc','ds1');
      filter = field.filters['default'];
      field._range = [2,6];
      filter.range([3,5]);
      done();
    })
    it("should change the upper range", (done) => {
      filter.rangeHigh(6);
      filter._range.should.be.a.Array();
      filter._range.should.match([3,6]);
      done();
    });
    it("should be bounded by the limits", (done) => {
      filter.rangeHigh(7);
      filter._range.should.be.a.Array();
      filter._range.should.match([3,6]);
      done();
    });
    it("should adjust the lower range if out of bounds", (done) => {
      filter.rangeHigh(2);
      filter._range.should.be.a.Array();
      filter._range.should.match([2,2]);
      done();
    });
    it("should return a number if called with no value", (done) => {
      let result = filter.rangeHigh();
      result.should.be.a.Number();
      result.should.equal(5);
      done();
    });
  });
  describe("rangeLow(value)", (done) => {
    let field;
    let filter;
    beforeEach((done) => {
      field = new Field('gc','ds1');
      filter = field.filters['default'];
      field._range = [2,6];
      filter.range([3,5]);
      done();
    })
    it("should change the lower range", (done) => {
      filter.rangeLow(4);
      filter._range.should.be.a.Array();
      filter._range.should.match([4,5]);
      done();
    });
    it("should be bounded by the limits", (done) => {
      filter.rangeLow(0);
      filter._range.should.be.a.Array();
      filter._range.should.match([2,5]);
      done();
    });
    it("should adjust the upper range if out of bounds", (done) => {
      filter.rangeLow(6);
      filter._range.should.be.a.Array();
      filter._range.should.match([6,6]);
      done();
    });
    it("should return a number if called with no value", (done) => {
      let result = filter.rangeLow();
      result.should.be.a.Number();
      result.should.equal(3);
      done();
    });
  });
  describe("rangePercent(arr)", (done) => {
    let field;
    let filter;
    beforeEach((done) => {
      field = new Field('gc','ds1',{range:[0,200]});
      filter = field.filters['default'];
      done();
    })
    it("should show range as a percentage if called without values", (done) => {
      filter.range([20,160]);
      let result = filter.rangePercent();
      result.should.be.a.Array();
      result.should.match([10,80]);
      done();
    });
    it("should support log scaled variables", (done) => {
      field.range([1,10000]);
      field.scaleType('scaleLog');
      filter.range([10,1000]);
      let result = filter.rangePercent();
      result.should.be.a.Array();
      result.should.match([25,75]);
      done();
    });
    it("should update the range", (done) => {
      filter.rangePercent([50,75]);
      let result = filter.range();
      result.should.be.a.Array();
      result.should.match([100,150]);
      done();
    });
  });
  describe("rangeHighPercent(value)", (done) => {
    let field;
    let filter;
    beforeEach((done) => {
      field = new Field('gc','ds1',{range:[0,200]});
      filter = field.filters['default'];
      done();
    })
    it("should show upper range as a percentage if called without value", (done) => {
      filter.range([20,160]);
      let result = filter.rangeHighPercent();
      result.should.be.a.Number();
      result.should.equal(80);
      done();
    });
    it("should update the upper range", (done) => {
      filter.rangeHighPercent(75);
      let result = filter.rangeHigh();
      result.should.be.a.Number();
      result.should.equal(150);
      done();
    });
  });
  describe("rangeLowPercent(value)", (done) => {
    let field;
    let filter;
    beforeEach((done) => {
      field = new Field('gc','ds1',{range:[0,200]});
      filter = field.filters['default'];
      done();
    })
    it("should show lower range as a percentage if called without value", (done) => {
      filter.range([20,160]);
      let result = filter.rangeLowPercent();
      result.should.be.a.Number();
      result.should.equal(10);
      done();
    });
    it("should update the lower range", (done) => {
      filter.rangeLowPercent(50);
      let result = filter.rangeLow();
      result.should.be.a.Number();
      result.should.equal(100);
      done();
    });
  });
  describe("active(bool)", (done) => {
    let filter;
    beforeEach((done) => {
      filter = new Filter('gc','ds1');
      done();
    })
    it("should return state of filter", (done) => {
      filter._active = true;
      let result = filter.active();
      result.should.be.True();
      filter._active = false;
      result = filter.active();
      result.should.be.False();
      done();
    });
    it("should update the filter state", (done) => {
      filter._active = true;
      let result = filter.active(false);
      result.should.be.False();
      done();
    });
  });
  describe("inclusive(bool)", (done) => {
    let filter;
    beforeEach((done) => {
      filter = new Filter('gc','ds1');
      done();
    })
    it("should return state of filter", (done) => {
      filter._inclusive = true;
      let result = filter.inclusive();
      result.should.be.True();
      filter._inclusive = false;
      result = filter.inclusive();
      result.should.be.False();
      done();
    });
    it("should update the filter state", (done) => {
      filter._inclusive = true;
      let result = filter.inclusive(false);
      result.should.be.False();
      done();
    });
  });
});
