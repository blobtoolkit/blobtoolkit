const should = require('should');
const request = require('supertest');
const express = require('express');
const app = require('../app');

describe('Routes:', () => {
  describe('dataset:', () => {
    describe('/api/v1/dataset/all', () => {
      it('should return well formatted json', (done) => {
        request(app)
          .get('/api/v1/dataset/all')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(200);
            res.header['content-type'].should.match(/application\/json/)
            res.body.should.be.a.Array();
            res.body[0].should.have.properties('id','name','description');
            res.body[0].id.should.be.a.String();
            done();
          });
      });
    });
    describe('/api/v1/dataset/id/:dataset_id', () => {
      it('should return well formatted json', (done) => {
        request(app)
          .get('/api/v1/dataset/id/ds1')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(200);
            res.header['content-type'].should.match(/application\/json/)
            res.body.should.be.a.Object();
            res.body.should.have.properties('id','name','description','records');
            res.body.records.should.be.a.Number();
            done();
          });
      });
      it('should return 404 for an invalid request', (done) => {
        request(app)
          .get('/api/v1/dataset/id/noID')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(404);
            done();
          });
      });
    });
    describe('/api/v1/dataset/id/:dataset_id/:key', () => {
      it('should return values', (done) => {
        request(app)
          .get('/api/v1/dataset/id/ds1/records')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(200);
            res.header['content-type'].should.match(/application\/json/)
            res.body.should.be.a.Number();
            res.body.should.be.equal(10);
            done();
          });
      });
      it('should return fields by id', (done) => {
        request(app)
          .get('/api/v1/dataset/id/ds1/gc')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(200);
            res.header['content-type'].should.match(/application\/json/)
            res.body.should.be.a.Object();
            res.body.should.have.properties('_range');
            res.body._range.should.be.a.Array();
            res.body._range[0].should.be.a.Number();
            done();
          });
      });
      it('should look for nested fields', (done) => {
        request(app)
          .get('/api/v1/dataset/id/ds1/bam0_cov')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(200);
            res.header['content-type'].should.match(/application\/json/)
            res.body.should.be.a.Object();
            res.body.should.have.properties('_name');
            res.body._name.should.be.a.String();
            res.body._name.should.be.match('bam0');
            done();
          });
      });
      it('should return 404 for an invalid request', (done) => {
        request(app)
          .get('/api/v1/dataset/id/ds1/absent')
          .end((err,res) => {
            if (err) return done(err);
            res.status.should.equal(404);
            done();
          });
      });
    });
  })
});
