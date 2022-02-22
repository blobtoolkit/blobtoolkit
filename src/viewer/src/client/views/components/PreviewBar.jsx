import React from 'react';

const PreviewBar = ({ css, x, y, width, height }) => (
  <g className={css} transform={'translate('+x+','+y+')'}>
    <rect x='1' width={width} height={height}/>
  </g>
);

export default PreviewBar;
