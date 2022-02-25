import React from 'react';
import { connect } from 'react-redux'
import { getTransformLines }  from '../reducers/plotOverlay'
import styles from './Plot.scss'

export default class PlotTransformLines extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        getTransformLines(state)
      )
    }
  }

  render(){
    const ConnectedTransformLines = connect(
      this.mapStateToProps
    )(TransformLines)
    return (
      <ConnectedTransformLines {...this.props}/>
    )
  }
}

const TransformLines = ({ lines = [], params, origins}) => {
  if (params.factor != 0){
    let origin
    if (params.origin.x == 0){
      origin = <g>
        {lines.map((line,i) =>
          <path d={line} key={i} className={styles.transform_line}  clipPath="url(#plot-area)"/>
        )}
      </g>
    }
    else {
      origin = <g>
        <path d={origins[1]} key={'x_origin'} className={styles.bold_transform_line}  clipPath="url(#plot-area)"/>
        <path d={origins[0]} key={'y_origin'} className={styles.bold_transform_line}  clipPath="url(#plot-area)"/>
      </g>
    }
    return (
      origin
    )
  }
  return null
}
