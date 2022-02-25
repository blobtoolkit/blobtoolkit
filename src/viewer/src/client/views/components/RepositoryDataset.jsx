import React from 'react'
import { connect } from 'react-redux'
import PropTypes from 'prop-types'
import styles from './Repository.scss';
import { Link } from 'react-router-dom'
import Spinner from './Spinner'
import { createSelector } from 'reselect'
import { getDatasetMeta } from '../reducers/repository'

class RepositoryDataset extends React.Component {
  constructor(props) {
    super(props);
  }

  componentDidMount(){
    this.props.onMount(this.props.id);
  }

  render(){
    const mapStateToProps = state => {
      return {
        meta: getDatasetMeta(state,this.props.id)
      }
    }
    const ConnectedDataset = connect(
      mapStateToProps,
      false
    )(Dataset)

    return (
      <Link to={'/view/'+this.props.id}>
        <div className={styles.dataset} onClick={()=>{this.props.onMount(this.props.id)}}>
          <ConnectedDataset {...this.props}/>
        </div>
      </Link>
    )
  }


}

class Dataset extends React.Component {
  render(){
    if (!this.props.meta){
      return <Spinner/>
    }
    return (
      <div>
        {this.props.meta.name}<br/>
        {this.props.meta.records} {this.props.meta.record_type}
      </div>
    )
  }
}

export default RepositoryDataset
