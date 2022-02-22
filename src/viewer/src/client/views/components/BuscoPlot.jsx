import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
import { getBuscoSets, getAllBuscoCSV } from '../reducers/summary'
import Busco from './Busco'
import BuscoData from './BuscoData'
import PlotMissing from './PlotMissing'
import getSelectedDatasetMeta from '../reducers/dataset'

class BuscoSets extends React.Component {

  shouldComponentUpdate(nextProps, nextState){
    if (Object.keys(this.props.buscoData).length === 0){
      return false
    }
    if (nextProps.buscoData.every((x,i)=>(this.props.buscoData[i] && x.length == this.props.buscoData[i].length),this)){
      return false
    }
    return false
  }

  render(){
    if (!this.props.buscoSets) {
      return (
        <PlotMissing view='busco' name={this.props.datasetName}/>
      )
    }
    let buscos = []
    this.props.buscoSets.forEach((id,i)=>{
      buscos.push(<Busco key={i} id={id}/>)
    })
    return (
      <div className={styles.busco_fill_parent}>
        <BuscoData {...this.props}/>
        {buscos}
      </div>
    )
  }
}

class BuscoPlot extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        buscoSets: getBuscoSets(state),
        buscoData: getAllBuscoCSV(state)
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

export default BuscoPlot
