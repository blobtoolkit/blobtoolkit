import React from 'react';
import styles from './Layout.scss'
import uoeLogo from './img/uoe-logo.png'
import ebiLogo from './img/ebi-logo.png'
import sangerLogo from './img/sanger-logo.png'
import bbsrcLogo from './img/bbsrc-logo.png'

const BTKLogos = ({}) => {
  return (
    <span className={styles.logo}>
      <a href='https://www.sanger.ac.uk'>
        <img src={sangerLogo} alt='Sanger' />
      </a>
      <a href='https://www.ebi.ac.uk'>
        <img src={ebiLogo} alt='EBI' />
      </a>
      <a href='https://www.ed.ac.uk'>
        <img src={uoeLogo} alt='UoE' />
      </a>
      <a href='https://bbsrc.ukri.org'>
        <img src={bbsrcLogo} alt='BBSRC' />
      </a>
    </span>
  )
}

export default BTKLogos
