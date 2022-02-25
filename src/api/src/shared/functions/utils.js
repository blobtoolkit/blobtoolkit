var Promise = require('promise');

const waitOn = async (value,promise) => {
  return await Promise.resolve(value);
}

// TODO: write test
const removeNestedKeys = (arr,keys,nestarr) => {
  arr.forEach((obj) => {
    nestarr.forEach((name) => {
      if (obj.hasOwnProperty(name)){
        removeNestedKeys(obj[name],keys,nestarr)
      }
    })
    keys.forEach((key) => {
      if (obj.hasOwnProperty(key)){
        delete obj[key];
      }
    });
  })
  return true;
}

const asKeyValue = (arr) => {
  let obj = {};
  let keys = [];
  let values = arr.map(function(v) {
    if (!obj.hasOwnProperty(v)) {
      obj[v] = keys.length;
      keys.push(v);
    }
    return obj[v];
  });
  return {values:values,keys:keys};
};

const groupValuesBy = (arr, by, val) => {
  grouped = arr.reduce((rv, x) => {
    (rv[x[by]] = rv[x[by]] || []).push(x[val]);
    return rv;
  }, {});
  return Promise.resolve(grouped)
};

const valueAtPath = (obj,path) => {
  path.forEach((p) => {
    obj = obj[p];
  })
  //return Promise.resolve(obj);
  return obj;
}

const newField = (id,options = {}) => {
  let field = {};
  field.id = id;
  options.name = options.name || id;
  Object.keys(options).forEach(function(key){
    field[key] = options[key];
  })
  return field;
}

const entryByKeyValue = (arr,key,value) => {
  let result = [];
  arr.forEach((obj)=>{
    if (obj.hasOwnProperty(key)){
      if (obj[key] == value){
        result.push(obj);
      }
    }
  })
  /*if (result.length == 1){
    result = result[0];
  }
  else */if (result.length == 0) {
    result = undefined;
  }
  return result;
}

const nestedEntryByKeyValue = (arr,key,value,nestarr) => {
  let result = entryByKeyValue(arr,key,value);
  if (!result){
    result = [];
    arr.forEach((obj) => {
      nestarr.forEach((name) => {
        if (obj.hasOwnProperty(name)){
          let res = nestedEntryByKeyValue(obj[name],key,value,nestarr)
          if (res){
            result = result.concat(res)
          }
        }
      })
    })
  }
  return result;
}

const objectToArray = (obj) => {
  let i=0;
  let arr=[];
  for (let o in obj){
    arr[i++] = obj[o];
  }
  return arr;
}


module.exports = {
  waitOn,
  asKeyValue,
  groupValuesBy,
  valueAtPath,
  newField,
  entryByKeyValue,
  nestedEntryByKeyValue,
  removeNestedKeys,
  objectToArray
}
