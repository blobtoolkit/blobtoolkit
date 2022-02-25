import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'

class StaticWarning extends React.Component {
  constructor(props) {
    super(props);
    this.state = { dismissed: false };
  }

  toggleDismissed(dismissed){
    this.setState({dismissed})
  }

  render(){
    if (this.state.dismissed){
      return <span className={styles.info} onClick={()=>this.toggleDismissed(false)}/>
    }
    else {
      let message
      if (this.props.threshold < this.props.records){
        message = 'The number of records in this dataset ('+(this.props.records*1).toLocaleString()+') is greater than the current static threshold ('+(this.props.threshold*1).toLocaleString()+').'
      }
      else {
        message = 'The view mode has been set to static.'
      }
      return (
        <div className={styles.warning}>
          <h2>This is a static view of dataset {this.props.name}</h2>
          <p>{message}</p>
          Use the Settings menu to switch to interactive mode.
          <span className={styles.dismiss} onClick={()=>this.toggleDismissed(true)}/>
        </div>
      )
    }

  }
}


export default StaticWarning
