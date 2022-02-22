import React from 'react'
import styles from './Filters.scss'
import { connect } from 'react-redux'
import { getFilteredBarsForFieldId } from '../reducers/preview'
import Spinner from './Spinner'
import PreviewBars from './PreviewBars'

class FilterPreview extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        bars: getFilteredBarsForFieldId(state, this.props.fieldId),
        barcss: styles.bar
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
      }
    }
  }

  render(){
    const ConnectedPreview = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(FilteredDataPreview)
    return <ConnectedPreview {...this.props}/>
  }
}

class FilteredDataPreview extends React.Component {

  render() {
    return (
      <g className={styles.filter_preview_container}>
        <PreviewBars bars={this.props.bars} barcss={this.props.barcss} />
      </g>
    );
  }

}

export default FilterPreview
