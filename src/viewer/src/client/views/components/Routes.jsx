import React from 'react'
import { connect } from 'react-redux'
import { Router, Switch, Route } from 'react-router-dom'
import { getDatasetIsActive } from '../reducers/repository'
import Dataset from './Dataset';
import Layout from './Layout';
import Repository from './Repository';
import { NotFound } from './Error';
import history from '../reducers/history';

export const Routes = (active) => {
  return (
    <Layout {...active}/>
  )
  // return (
  //   <Router history={history}>
  //     <Switch>
  //       <Route path="/dataset/:datasetId/:view?" render={(props)=>(<Layout {...props}/>)}/>
  //       <Route path="/:searchTerm/dataset/:datasetId/:view?" render={(props)=>(<Layout {...props}/>)}/>
  //       <Route path="/:searchTerm?" render={(props)=>(<Layout {...props}/>)}/>
  //     </Switch>
  //   </Router>
  // )
}

export default Routes
