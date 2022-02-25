import React from 'react'
import styles from './Menu.scss';

class MenuItem extends React.Component {
  render(){
    return (
      <div className={styles.outer}>
        <div className={styles.header}>
          <h1 className={styles.inline}>{this.props.name}</h1>
        </div>
      </div>
    )
  }
}


export default MenuItem
