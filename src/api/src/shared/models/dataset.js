const Field = require('./field');
const utils = require('../functions/utils');

function Dataset(id) {
  this.id = id;
  this.fields = {};
  this.hierarchy = {};
};

module.exports = Dataset;

const addFields = function(fields,m = {},parent) {
  fields.forEach((f) => {
    let id = f.id || f._id;
    let p;
    if (!parent){
      this.hierarchy[id] = {};
      p = this.hierarchy[id];
    }
    else {
      parent[id] = {};
      p = parent[id];
    }
    let meta = {};
    Object.keys(m).forEach((key) => {
      if (key != 'children' && key != 'data'){
        meta[key] = m[key];
      }
    })
    Object.keys(f).forEach((key) => {
      if (key != 'children' && key != 'data'){
        meta[key] = f[key];
      }
    })
    if (f.children || f._children){
      this.addFields(f.children || f._children,meta,p);
    }
    else {
      if (f.preload || f._preload){
        f._active = f.preload || f._preload;
      }
      this.addField(f.id || f._id,meta);
      if (f.data || f._data){
        this.addFields(f.data || f._data,meta,p);
      }
    }
  })
}

const addField = function(fId,meta){
  if (!this.fields.hasOwnProperty(fId)){
    this.fields[fId] = new Field(fId,this,meta);
  }
  return this.fields[fId];
}

Dataset.prototype = {
  addFields,
  addField
}
