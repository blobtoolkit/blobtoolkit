import React from 'react'
import styles from './Plot.scss'
import PlotSideBinsSVG from './PlotSideBinsSVG'

export default class XPlot extends React.Component {
  render(){
    if (this.props.plotGraphics == 'canvas'){
        //plotContainer = <PlotSquareBinsCanvas />
    }
    else {
      return (
        <div className={styles.x_axis_plot}>
          <svg ref={(elem) => { this.svg = elem; }}
            className={styles.axis_plot}
            viewBox='0 0 1000 300'>
            <PlotSideBinsSVG axis='x'/>
          </svg>
        </div>

      )
    }
  }
}
