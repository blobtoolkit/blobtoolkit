import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'

export class NoHitWarning extends React.Component {
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
    else if (this.props.circleLimit){
      return (
        <div className={styles.warning}>
          <h2>'no-hit' data are not shown</h2>
          This dataset has over {Number(this.props.circleLimit).toLocaleString()} records so 'no-hit' circles are not plotted.
          <span className={styles.dismiss} onClick={()=>this.toggleDismissed(true)}/>
        </div>
      )
    }
    else if (this.props.nohitThreshold){
      return (
        <div className={styles.warning}>
          <h2>'no-hit' data are not shown</h2>
          This dataset has over {Number(this.props.nohitThreshold).toLocaleString()} records so "no-hit" records have been filtered out.
          Use the Filters menu to restore the hidden category data.
          <span className={styles.dismiss} onClick={()=>this.toggleDismissed(true)}/>
        </div>
      )
    }
    return null

  }
}


export default NoHitWarning
