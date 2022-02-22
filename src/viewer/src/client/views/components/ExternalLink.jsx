import React from 'react';
import styles from './Layout.scss'
import { Link } from 'react-router'

const ExternalLink = ({ title, url }) => {
  return (
    <a
      className={styles.external_link}
      href={url}
      target={title}
    >
    {title}
    </a>
  )
}

export default ExternalLink
