import React from 'react'
import { connect } from 'react-redux'
import { getCategoryListForMainPlot,
  getScatterPlotDataForCategoryIndex,
  getSquareBinPlotDataForCategoryIndex,
  getHexBinPlotDataForCategoryIndex } from '../reducers/plotData'
import PlotBubblesCanvas from './PlotBubblesCanvas'
import PlotSquareBinsSVG from './PlotSquareBinsSVG'
import PlotHexBinsSVG from './PlotHexBinsSVG'
import styles from './Plot.scss'

class PlotLayer extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        //return getScatterPlotDataForCategoryIndex(state,props.index)
        //return getScatterCanvasForCategoryIndex(state,props.index)
        //return getSquareBinPlotDataForCategoryIndex(state,props.index)
        return getHexBinPlotDataForCategoryIndex(state,props.index)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
      }
    }
  }

  render(){
    const ConnectedPlotLayer = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(CurrentLayer)
    return (
      <ConnectedPlotLayer {...this.props}/>
    )
  }
}

// const CurrentLayer = ({ data, bins, color }) => (
//   <a>{bins.id}</a>
// );

class CurrentLayer extends React.Component {
  render(){
    if (this.props.type == 'bubblesCanvas'){
      return (
        <div className={styles.fill_parent} style={{zIndex:this.props.zIndex}}>
          <PlotBubblesCanvas index={this.props.index} bubbles={this.props.data || []} color={this.props.color} />
        </div>
      )
    }
    else if (this.props.type == 'squareBinsSVG'){
      return (
        <g className={styles.fill_parent} style={{zIndex:this.props.zIndex}}>
          <PlotSquareBinsSVG index={this.props.index} data={this.props.data || []} color={this.props.color} side={this.props.side || 1} />
        </g>
      )
    }
    else if (this.props.type == 'hexBinsSVG'){
      return (
        <g className={styles.fill_parent} style={{zIndex:this.props.zIndex}}>
          <PlotHexBinsSVG index={this.props.index} data={this.props.data || []} color={this.props.color} />
        </g>
      )
    }
  }
}

export default PlotLayer
