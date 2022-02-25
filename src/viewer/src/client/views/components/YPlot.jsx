import React from 'react'
import styles from './Plot.scss'
import PlotSideBinsSVG from './PlotSideBinsSVG'

export default class YPlot extends React.Component {
  render(){
    if (this.props.plotGraphics == 'canvas'){
        //plotContainer = <PlotSquareBinsCanvas />
    }
    else {
      return (
        <div className={styles.y_axis_plot}>
          <svg ref={(elem) => { this.svg = elem; }}
            className={styles.axis_plot}
            viewBox='0 0 1000 300'>
            <PlotSideBinsSVG axis='y'/>
          </svg>
        </div>

      )
    }
  }
}
