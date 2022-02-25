import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'

class PlotMissing extends React.Component {
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
          <h2>This view is not yet available for dataset {this.props.name}</h2>
          <p>No data are available for the {this.props.view} view as we have not yet run the analysis.</p>
          Use the Settings menu to choose an alternative view.
          <span className={styles.dismiss} onClick={()=>this.toggleDismissed(true)}/>
        </div>
      )
    }

  }
}


export default PlotMissing
