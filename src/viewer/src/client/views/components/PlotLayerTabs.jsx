import React from 'react';
import styles from './Plot.scss'

const PlotLayerTabs = ({ children }) => (
    <div className={styles.layer_tab_holder}>
      {children}
    </div>
);

export default PlotLayerTabs;
