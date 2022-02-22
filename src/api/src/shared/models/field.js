const d3 = require('d3')
const Filter = require('./filter');

function Field(id,dataset,meta) {
  this._id = id;
  this.dataset = dataset;
  if (meta){
    Object.keys(meta).forEach((key)=>{
      if (key == 'range'){
        this.range(meta[key]);
      }
      else {
        let newkey = key.match(/^_/) ? key : '_'+key
        this[newkey] = meta[key];
      }
    })
  }
  if (this._range){
    this.scaleType(this._scale || 'scaleLinear');
  }
  let filter = new Filter('default',this);
  this.filters = {default:filter};
  if (this._range){
    filter._limits = this._range;
    filter.scale = this.scale;
  }
};

module.exports = Field;

const allowedScales = {
    scaleLinear: 'scaleLinear',
    scaleLog: 'scaleLog',
    scaleSqrt: 'scaleSqrt',
    scalePow: 'scalePow',
    scaleOrdinal: 'scaleOrdinal'
};

const scaleType = function(value){
  let result
  if (typeof value === 'string'){
    if (allowedScales[value]){
      let domain = this.scale ? this.scale.domain() : this._range;
      this.scale = d3[allowedScales[value]]().domain(domain).range([0,100]);
      this.scale.scaleType = allowedScales[value];
      result = this.scale.scaleType;
    }
    else {
      result = undefined;
    }
  }
  else {
    result = this.scale.scaleType;
  }
  return result;
}

const range = function(arr){
  if (Array.isArray(arr) && arr.length > 0){
    this._range = [Math.min(...arr),Math.max(...arr)];
    if (this.scale){
      this.scale.domain(this._range);
    }
  }
  return this._range;
}

const rangeHigh = function(value){
  if (typeof value === 'number') {
    this._range[1] = value;
    this._range[0] = Math.min(this._range[0],value);
    this.scale.domain(this._range);
  }
  return this._range[1];
}

const rangeLow = function(value){
  if (typeof value === 'number') {
    this._range[0] = value;
    this._range[1] = Math.max(this._range[1],value);
    this.scale.domain(this._range);
  }
  return this._range[0];
}

const filterToList = async function(arr){
  let low = this.filters['default'].rangeLow();
  let high = this.filters['default'].rangeHigh();
  let inclusive = this.filters['default'].inclusive();
  let filter;
  if (inclusive){
    filter = (value) => { return value >= low && value <= high }
  }
  else {
    filter = (value) => { return value <= low || value >= high }
  }
  let result = await this.loadData().then(()=>{
    let data = this.data;
    let values = data.values;
    let ret = [];
    if (arr){
      arr.forEach(function(i){
        if (filter(values[i])){
          ret.push(i);
        }
      })
    }
    else {
      values.forEach(function(v,i){
        if (filter(v)){
          ret.push(i);
        }
      })
    }
    return ret;
  })
  return Promise.resolve(result)
}

const applyListFilter = async function(arr){
  let result = await this.loadData().then(()=>{
    let data = this.data;
    let keys = data.keys;
    let values = data.values;
    let ret = [];
    arr.forEach(function(i){
      if (keys){
        ret.push(keys[values[i]]);
      }
      else {
        ret.push(values[i]);
      }
    })
    this.filtered = ret;
    return ret;
  })
  return Promise.resolve(result)
}

Field.prototype = {
  scaleType,
  range,
  rangeHigh,
  rangeLow,
  filterToList,
  applyListFilter
}
