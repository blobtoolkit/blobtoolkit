import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
import { CSVLink, CSVDownload } from 'react-csv'
import { getDatasetCSVdata } from '../reducers/datasetTable'
import { getSearchTerm } from '../reducers/location'

export const DatasetCSVComponent = ({data,term}) => {
  if (!data || data.length == 0) return null
  return (
    <CSVDownload
        filename={term+'.csv'}
        data={data}
        target='_blank'/>
  )
}


class DatasetTableCSV extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        data: getDatasetCSVdata(state),
        term: getSearchTerm(state)
      }
    }
  }

  render(){
    const ConnectedCSV = connect(
      this.mapStateToProps
    )(DatasetCSVComponent)
    return <ConnectedCSV {...this.props}/>
  }
}

export default DatasetTableCSV
