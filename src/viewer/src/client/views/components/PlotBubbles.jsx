import React from 'react';
import PlotBubble from './PlotBubble'

const PlotBubbles = ({ bubbles, bubblecss }) => (
    <g>
      {bubbles.map(bubble =>
        <PlotBubble key={bubble.id} {...bubble} css={bubblecss} />
      )}
    </g>
);

export default PlotBubbles;
