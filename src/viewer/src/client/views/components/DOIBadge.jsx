import React from 'react';
import styles from './Layout.scss'

const DOIBadge = ({}) => {
  return (
    <div className={styles.doi_badge}>
      <a href='https://doi.org/10.5281/zenodo.1134794'>
        <img src='https://img.shields.io/badge/DOI-10.5281%2Fzenodo.1134794-9c528b?style=flat&labelColor=414a51' alt='DOI' />
      </a>
    </div>

  )
}

export default DOIBadge
