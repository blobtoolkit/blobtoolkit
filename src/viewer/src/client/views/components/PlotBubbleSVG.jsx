import React from 'react';
import { connect } from 'react-redux'

const PlotBubbleSVG = ({ x, y, r, plotShape }) => (
  <circle cx={x} cy={y} r={r} style={plotShape.circle}/>
);

export default PlotBubbleSVG
