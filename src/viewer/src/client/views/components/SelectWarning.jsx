import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'
import { getActiveSelection } from '../reducers/filter'
import { getRecordCount } from '../reducers/summary'


class Warning extends React.Component {
  constructor(props) {
    super(props);
    this.state = { dismissed: true };
  }

  toggleDismissed(dismissed){
    this.setState({dismissed})
  }

  render(){
    let selCount = this.props.activeSelection ? this.props.activeSelection.length : 0
    if (selCount == 0 || selCount == this.props.recordCount) return null
    if (this.state.dismissed){
      return <span className={styles.info} onClick={()=>this.toggleDismissed(false)}/>
    }
    else {
      return (
        <div className={styles.warning}>
          <h2>Selection-based filtering is not captured in the URL</h2>
          <p>This dataset is filtered by an active selection that cannot be captured in the URL.</p>
          Export the current settings as a list to allow this session to be reproduced.
          <span className={styles.dismiss} onClick={()=>this.toggleDismissed(true)}/>
        </div>
      )
    }

  }
}

class SelectWarning extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        activeSelection: getActiveSelection(state),
        recordCount: getRecordCount(state)
      }
    }
  }

  render(){
    const ConnectedWarning = connect(
      this.mapStateToProps,
    )(Warning)
    return <ConnectedWarning {...this.props}/>
  }
}

export default SelectWarning
