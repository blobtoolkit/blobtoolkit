import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss';
import Spinner from './Spinner'

class MenuDisplaySet extends React.Component {
  render(){
    return (
      <div>
        {this.props.children}
      </div>
    )
  }
}


export default MenuDisplaySet
