import React from 'react'
import { connect } from 'react-redux'
import {
  filterToList
} from '../reducers/filter'
import {
  getFilteredSummary
} from '../reducers/preview'
import Spinner from './Spinner'
import styles from './Datasets.scss'


const mapStateToProps = state => getFilteredSummary(state)

const mapDispatchToProps = dispatch => {
  return {
    onClick: () => dispatch(filterToList())
  }
}

class ApplyFilters extends React.Component {

  render(){
    let text = this.props.selected+'/'+this.props.count+' ('+this.props.percentage+') '+this.props.type + 's'
    //let percent = this.props.selected/this.props.count*100
    return (
      <div className={styles.dataset_controls}>
        <div className={styles.detail_container}
          data-tip data-place='bottom' data-for='filter-summary'>
          {text}
          <div className={styles.detail} style={{width:this.props.percentage}}>
            {text}
          </div>
        </div>

      </div>
    )
  }
}

const DatasetApplyFilters = connect(
  mapStateToProps,
  mapDispatchToProps
)(ApplyFilters)

export default DatasetApplyFilters
