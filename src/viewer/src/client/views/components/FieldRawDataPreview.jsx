import React from 'react'
import styles from './Fields.scss'
import { connect } from 'react-redux'
import { getBarsForFieldId } from '../reducers/field'
import { setDimension } from '../reducers/dimension'
import Spinner from './Spinner'
import PreviewBars from './PreviewBars'

class FieldRawDataPreview extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        bars: getBarsForFieldId(state, this.props.fieldId),
        barcss: styles.bar
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        onMount: (obj) => {}//dispatch(setDimension({id:"preview",height:300,width:1225}))
      }
    }
  }

  render(){
    let previewType = RangeDataPreview;
    if (this.props.fieldId.match(/[ms]$/)){
    //  previewType = CatDataPreview;
    }
    const ConnectedRawDataPreview = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(previewType)
    return <ConnectedRawDataPreview {...this.props}/>
  }
}

class RangeDataPreview extends React.Component {

  componentDidMount(){
    this.props.onMount({
      id:'preview',
      width:1000,
      height:250
    })
  }

  render() {
    let divider
    if (this.props.meta.clamp && this.props.meta.limit[0] < this.props.meta.clamp){
      divider = <line className={styles.clamped_divider}
                      x1={25} x2={25} y1={0} y2={106.67}/>
    }
    return (
      <g className={styles.data_preview_container}>
        {divider}
        <PreviewBars bars={this.props.bars} barcss={this.props.barcss} />
      </g>
    );
  }

}

class CatDataPreview extends React.Component {

  componentDidMount(){
    this.props.onMount({
      id:'preview',
      width:this.svg.clientWidth,
      height:this.svg.clientHeight
    })
  }

  render() {
    return (
      <div className={styles.category_data_preview_container}>
        <svg ref={(elem) => { this.svg = elem; }}>
          <PreviewBars bars={this.props.bars} barcss={this.props.barcss} />
        </svg>
      </div>
    );
  }

}


export default FieldRawDataPreview
