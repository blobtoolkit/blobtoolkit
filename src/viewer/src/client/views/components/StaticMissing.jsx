import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'

class StaticMissing extends React.Component {
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
      return (
        <div className={styles.warning}>
          <h2>This view is not available for dataset {this.props.name}</h2>
          The {this.props.view} view cannot be loaded in static mode, choose a different view or use the Settings menu to switch to interactive mode.
          <span className={styles.dismiss} onClick={()=>this.toggleDismissed(true)}/>
        </div>
      )
    }

  }
}


export default StaticMissing
