function Repository(id){
  this.id = id;
};

module.exports = Repository;

const getId = function() {
  return this.id || 'nothing';
}

Repository.prototype = {
  getId
}
