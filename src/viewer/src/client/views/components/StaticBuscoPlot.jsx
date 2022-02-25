import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
// import { getBuscoSets, getAllBuscoCSV } from '../reducers/summary'
import StaticBusco from './StaticBusco'
import StaticBuscoData from './StaticBuscoData'
import PlotMissing from './PlotMissing'
import { getDatasetID } from '../reducers/location'
import { getSelectedDatasetMeta } from '../reducers/dataset'

class BuscoSets extends React.Component {

  render(){
    let buscoSets = this.props.meta.summaryStats.busco
    if (!buscoSets) {
      return (
        <PlotMissing view='busco' name={this.props.datasetId}/>
      )
    }
    let buscos = []
    Object.keys(buscoSets).forEach((id,i)=>{
      buscos.push(<StaticBusco key={i} id={id} data={buscoSets[id]}/>)
    })
    // <BuscoData {...this.props}/>
    return (
      <div className={styles.busco_fill_parent} style={{overflow:'scroll'}}>
        <StaticBuscoData {...this.props}/>
        {buscos}
      </div>
    )
  }
}

class StaticBuscoPlot extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        meta: getSelectedDatasetMeta(state),
        datasetId: getDatasetID(state)
      }
    }
  }

  render(){
    const ConnectedBuscoSets = connect(
      this.mapStateToProps
    )(BuscoSets)
    return <ConnectedBuscoSets {...this.props}/>
  }
}

export default StaticBuscoPlot
