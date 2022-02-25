const d3 = require('d3')

function Filter(id,field,meta) {
  this._id = id;
  this.field = field;
  if (meta){
    Object.keys(meta).forEach((key)=>{
      this[key] = meta[key];
    })
  }
  this._name = this._name || this.id;
  this._active = false;
  this._inclusive = true;
  if (this.field._range){
    this.range(this.field._range);
  }
  //this.scaleType('scaleLinear');
};

module.exports = Filter;

const f3 = (value) => {
  return Number(d3.format('.3')(value));
}

const range = function(arr){
  if (Array.isArray(arr) && arr.length > 0){
    this._range = [Math.min(...arr),Math.max(...arr,)];
    this._range[0] = Math.max(this._range[0],this.field.rangeLow());
    this._range[1] = Math.min(this._range[1],this.field.rangeHigh());
  }
  return this._range;
}

const rangeHigh = function(value){
  if (typeof value === 'number') {
    this.range([Math.min(this._range[0],value),value]);
  }
  return this._range[1];
}

const rangeLow = function(value){
  if (typeof value === 'number') {
    this.range([value,Math.max(this._range[1],value)]);
  }
  return this._range[0];
}

const rangePercent = function(arr){
  let scale = this.field.scale;
  if (Array.isArray(arr) && arr.length > 0) {
    this.range(arr.map((a)=>{return scale.invert(a)}))
  }
  else {
    arr = [f3(scale(this._range[0])),f3(scale(this._range[1]))];
  }
  return arr;
}

const rangeHighPercent = function(value){
  let scale = this.field.scale;
  if (typeof value === 'number') {
    this.rangeHigh(scale.invert(value))
  }
  else {
    value = f3(scale(this._range[1]));
  }
  return value;
}

const rangeLowPercent = function(value){
  let scale = this.field.scale;
  if (typeof value === 'number') {
    this.rangeLow(scale.invert(value))
  }
  else {
    value = f3(scale(this._range[0]));
  }
  return value;
}

const active = function(bool){
  if (typeof bool != 'undefined') {
    if (bool){
      this._active = true;
    }
    else {
      this._active = false;
    }
  }
  return this._active;
}

const inclusive = function(bool){
  if (typeof bool != 'undefined') {
    if (bool){
      this._inclusive = true;
    }
    else {
      this._inclusive = false;
    }
  }
  return this._inclusive;
}


Filter.prototype = {
  range,
  rangeHigh,
  rangeLow,
  rangePercent,
  rangeHighPercent,
  rangeLowPercent,
  active,
  inclusive
}
