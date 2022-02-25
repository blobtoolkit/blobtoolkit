import React from 'react';

const PlotHexBinSVG = ({ x, y, zs, ids, points, color = '#999' }) => {
  //let width = Math.sqrt(zs.reduce((a,b)=>Math.max(a,b)))*2
  return (
    <polygon className='hex' fill={color} stroke={color} points={points}/>
  )
};

export default PlotHexBinSVG;
