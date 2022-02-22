import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'

export class NoBlobWarning extends React.Component {
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
          <h2>Blob view is not available for this dataset</h2>
          No coverage information is available for this dataset.
          {this.props.source && <span><br/> {this.props.source} view is shown instead.</span>}
          <p>Use the Settings menu to choose a different view.</p>
          <span className={styles.dismiss} onClick={()=>this.toggleDismissed(true)}/>
        </div>
      )
    }

  }
}


export default NoBlobWarning
