import React from 'react';
import styles from './Icon.scss'

const TextIcon = ({ title, invert, active, onIconClick = ()=>{} }) => {
  let css = styles.text_icon
  if (active) css += ' '+styles.active
  if (invert) css += ' '+styles.invert
  return (
    <a className={css} onClick={onIconClick}>{title}</a>
  )
}

export default TextIcon
