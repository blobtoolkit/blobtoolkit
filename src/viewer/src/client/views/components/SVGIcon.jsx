import React from 'react';
import styles from './Icon.scss'

const SVGIcon = ({ sprite, active, invert, onIconClick = ()=>{} }) => {
  let css = styles.text_icon
  if (active) css += ' '+styles.active
  if (invert) css += ' '+styles.invert
  return (
    <div className={styles.icon}
      data-tip data-for={sprite.id}>
      <svg viewBox={sprite.viewBox}
        className={css}
        onClick={onIconClick}>
        <use xlinkHref={'#'+sprite.id} />
      </svg>
    </div>
  )
}

export default SVGIcon
