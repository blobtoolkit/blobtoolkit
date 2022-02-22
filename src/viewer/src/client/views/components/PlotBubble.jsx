import React from 'react';

const PlotBubble = ({ css, cx, cy, r }) => (
  <circle cx={cx} cy={cy} r={r} className={css}/>
);

export default PlotBubble;
