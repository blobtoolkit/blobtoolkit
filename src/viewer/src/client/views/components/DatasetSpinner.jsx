import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'
import { loadDataset, getDatasetIsActive } from '../reducers/repository'
import { getView, getDatasetID } from '../reducers/location'
import GetStarted from './GetStarted'
import MainPlot from './MainPlot'
import CumulativePlot from './CumulativePlot'
import DetailPlot from './DetailPlot'
import SnailPlot from './SnailPlot'
import TablePlot from './TablePlot'
import TreeMapPlot from './TreeMapPlot'
import Spinner from './Spinner'

class DatasetSpinner extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        active: getDatasetIsActive(state),
        id: getDatasetID(state)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        onLoad: (id) => dispatch(loadDataset(id))
      }
    }
  }

  render(){
    const ConnectedSpinner = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(Spinner)
    return <ConnectedSpinner {...this.props}/>
  }
}

export default DatasetSpinner
