import React from 'react';
import { connect } from 'react-redux'
import styles from './Plot.scss'
import PlotBubbleSVG from './PlotBubbleSVG'
import { getCirclePlotDataForCategoryIndex }  from '../reducers/plotData'
import { getPlotResolution, getZScale }  from '../reducers/plotParameters'
import { plotShapes } from '../reducers/plotStyles'

export default class PlotBubblesSVG extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        let val = getCirclePlotDataForCategoryIndex(state,this.props.x0)
        let obj = {}
        obj.color = val.color
        obj.bubbles = val.data
        obj.scale = getZScale(state)
        obj.res = getPlotResolution(state)
        obj.plotShape = plotShapes(state)
        return obj
      }
    }
  }

  render(){
    const ConnectedBubblesSVG = connect(
      this.mapStateToProps
    )(BubblesSVG)
    return (
      <ConnectedBubblesSVG {...this.props}/>
    )
  }
}

const BubblesSVG = ({ bubbles, color, plotShape }) => (
  <g fill={color}>
    {bubbles.map(bubble =>
      <PlotBubbleSVG key={bubble.id} {...bubble} plotShape={plotShape} />
    )}
  </g>
);
