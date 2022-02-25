import React from 'react';
import { connect } from 'react-redux'
import styles from './Plot.scss'
import PlotBubblesCanvas from './PlotBubblesCanvas'
import { getCategoryListForMainPlot }  from '../reducers/plotData'

export default class PlotBubblesCanvasLayers extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        getCategoryListForMainPlot(state)
      )
    }
  }

  render(){
    const ConnectedBubbleCanvasLayers = connect(
      this.mapStateToProps
    )(BubblesCanvasLayers)
    return (
      <ConnectedBubbleCanvasLayers {...this.props}/>
    )
  }
}

const BubblesCanvasLayers = ({ bins, records, circleLimit }) => (
    <div className={styles.fill_parent}>
      {
        bins.map((layer,i) => {
          if (records > circleLimit && layer.id == 'no-hit'){
            return null
          }
          return <PlotBubblesCanvas key={layer.id} x0={layer.x0} />
        })
      }
    </div>
);
