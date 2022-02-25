import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
import { CSVLink, CSVDownload } from 'react-csv'
import { getCSVdata } from '../reducers/summary'
import { getDatasetID } from '../reducers/location'

export const DownloadCSVComponent = ({data,dataset}) => {
  if (!data || data.length == 0) return null
  return (
    <CSVDownload
        filename={dataset+'.csv'}
        data={data}
        target='_blank'/>
  )
}


class DownloadCSV extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        data: getCSVdata(state),
        dataset: getDatasetID(state)
      }
    }
  }

  render(){
    const ConnectedCSV = connect(
      this.mapStateToProps
    )(DownloadCSVComponent)
    return <ConnectedCSV {...this.props}/>
  }
}

export default DownloadCSV
