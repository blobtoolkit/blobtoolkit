import React from 'react';
import styles from './Plot.scss'

const PlotLayerTab = ({ layer, color, bin, onMouseOver }) => (
  <div className={styles.layer_tab} style={{backgroundColor:color}} onMouseEnter={onMouseOver}>
  <div className={styles.legend}>{bin.id}</div>
  </div>
);

export default PlotLayerTab;
