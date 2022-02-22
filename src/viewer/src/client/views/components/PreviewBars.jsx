import React from 'react';
import PreviewBar from './PreviewBar'

const PreviewBars = ({ bars, barcss }) => {
  if (bars) {
    return (
        <g>
          {bars.map(bar =>
            <PreviewBar key={bar.id} {...bar} css={barcss} />
          )}
        </g>
    )
  }
  return <g></g>
}


export default PreviewBars;
