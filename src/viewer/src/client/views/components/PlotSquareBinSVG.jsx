import React from 'react';

const PlotSquareBinSVG = ({ x, y, height, width, color = '#999' }) => {
  return (
    <rect fill={color} stroke={color} x={x} y={y} height={height} width={width} />
  )
};

export default PlotSquareBinSVG;
