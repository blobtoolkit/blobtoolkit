import React from 'react'
import { connect } from 'react-redux'
import { Router, Switch, Route } from 'react-router-dom'
import Dataset from './Dataset';
import Routes from './Routes';
import NoPlot from './NoPlot'
import Repository from './Repository';
import { NotFound } from './Error';
import history from '../reducers/history';
import { cookieConsent, analytics } from '../reducers/tracking';
import Spinner from './Spinner'
import { getDatasetID } from '../reducers/location'
import { getApiStatus, getDatasetIsActive } from '../reducers/repository'
import { getCookieConsent, getAnalytics, setCookieConsent, startAnalytics } from '../reducers/tracking'
import { fetchRepository, getRepositoryIsInitialised, getRepositoryIsFetching } from '../reducers/repository'
import ReactGA from 'react-ga';

class MainDiv extends React.Component {
  constructor(props) {
    super(props);
  }

  componentDidMount(){
    if (!this.props.initialised){
      this.props.onLoad()
    }
    if (!this.props.cookieConsent){
      if (this.props.cookies.cookies && String(this.props.cookies.cookies.CookieConsent) == 'true'){
        this.props.setCookieConsent(true)
      }
    }
    else {
      if (!this.props.analytics){
        this.props.startAnalytics(true)
      }
    }
  }

  componentDidUpdate(){
    if (this.props.cookieConsent){
      this.props.startAnalytics(true)
    }
  }

  render(){
    if (!this.props.apiStatus){
      console.log('error with API')
      return (<NoPlot reason='api'/>)
    }
    if (!this.props.initialised){
       return null
    }
    if (this.props.fetching){
      return <Spinner/>
    }
    return (
      <div>
        <Routes/>
      </div>
    )
  }
}


export default class Main extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => (
      {
        apiStatus: getApiStatus(state),
        initialised: getRepositoryIsInitialised(state),
        fetching: getRepositoryIsFetching(state),
        datasetId: getDatasetID(state),
        cookieConsent: getCookieConsent(state),
        analytics: getAnalytics(state)
      }
    )
    this.mapDispatchToProps = dispatch => (
      {
        onLoad: (searchTerm) => dispatch(fetchRepository(searchTerm)),
        setCookieConsent: (bool) => dispatch(setCookieConsent(bool)),
        startAnalytics: () => dispatch(startAnalytics())
      }
    )

  }

  render(){
    const ConnectedMain = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(MainDiv)
    return (
      <ConnectedMain {...this.props}/>
    )
  }
}

// module.exports = Main;
