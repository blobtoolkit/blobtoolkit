const should = require('should');
const request = require('supertest');
const express = require('express');
const app = require('../app');

describe('Routes:', () => {
  describe('field:', () => {
    describe('/api/v1/field/:dataset_id/:field_id', () => {
      it('should return well formatted json', (done) => {
        request(app)
          .get('/api/v1/field/ds1/gc')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(200);
            res.header['content-type'].should.match(/application\/json/)
            res.body.should.be.a.Object();
            res.body.should.have.property('values');
            res.body.values.should.be.a.Array();
            res.body.values.length.should.equal(10);
            res.body.values[0].should.be.a.Number();
            done();
          });
      });
    });
  });
});
