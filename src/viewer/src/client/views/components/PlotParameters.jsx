import React from 'react';
import { connect } from 'react-redux'
import styles from './Plot.scss'
import { getAxisTitle } from '../reducers/plotData'
import { getPlotShape,
  setPlotShape,
  getPlotResolution,
  setPlotResolution,
  getZReducer,
  setZReducer,
  getZScale,
  setZScale,
  getTransformFunctionParams,
  setTransformFunction } from '../reducers/plotParameters'

export default class PlotParameters extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        return {
          title:getAxisTitle(state,'z'),
          shape:getPlotShape(state),
          resolution:getPlotResolution(state),
          reducer:getZReducer(state),
          scale:getZScale(state),
          transform:getTransformFunctionParams(state)
        }
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        onSelectShape: shape => dispatch(setPlotShape(shape)),
        onChangeResolution: value => dispatch(setPlotResolution(value)),
        onSelectReducer: reducer => dispatch(setZReducer(reducer)),
        onSelectScale: reducer => dispatch(setZScale(reducer)),
        onChangeTransform: object => dispatch(setTransformFunction(object))
      }
    }
  }

  render(){
    const ConnectedPlotParameters = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(Parameters)
    return (
      <ConnectedPlotParameters {...this.props}/>
    )
  }
}
const Parameters = ({
  title,
  shape, onSelectShape,
  resolution, onChangeResolution,
  reducer, onSelectReducer,
  scale, onSelectScale,
  transform, onChangeTransform }) => {
  return (
    <div className={styles.axis_title+' '+styles.bottom_axis_title}>
      <span onClick={()=>onSelectShape('square')} className={shape == 'square' ? styles.active : ''}>square</span>
      <span onClick={()=>onSelectShape('hex')} className={shape == 'hex' ? styles.active : ''}>hex</span>
      <span onClick={()=>onSelectShape('circle')} className={shape == 'circle' ? styles.active : ''}>circle</span>
      <input onChange={(e)=>onChangeResolution(e.target.value)} type="range" value={resolution} min="5" max="50" step="1" className={styles.flip_horiz}/>
      <span onClick={()=>onSelectReducer('sum')} className={reducer.id == 'sum' ? styles.active : ''}>sum</span>
      <span onClick={()=>onSelectReducer('max')} className={reducer.id == 'max' ? styles.active : ''}>max</span>
      <span onClick={()=>onSelectReducer('min')} className={reducer.id == 'min' ? styles.active : ''}>min</span>
      <span onClick={()=>onSelectReducer('count')} className={reducer.id == 'count' ? styles.active : ''}>count</span>
      <span onClick={()=>onSelectReducer('mean')} className={reducer.id == 'mean' ? styles.active : ''}>mean</span>
      <span onClick={()=>onSelectScale('scaleLog')} className={scale == 'scaleLog' ? styles.active : ''}>log</span>
      <span onClick={()=>onSelectScale('scaleLinear')} className={scale == 'scaleLinear' ? styles.active : ''}>linear</span>
      <span onClick={()=>onSelectScale('scaleSqrt')} className={scale == 'scaleSqrt' ? styles.active : ''}>sqrt</span>
      <br/>
      <label htmlFor='transform_x'>x position: {transform.x} </label>
      <br/>
      <input id='transform_x' onChange={(e)=>onChangeTransform({x:e.target.value})} type="range" value={transform.x} min="0" max="1000" step="50"/>
      <br/>
      <label htmlFor='transform_order'>order: {transform.order} </label>
      <br/>
      <input id='transform_order' onChange={(e)=>onChangeTransform({order:e.target.value})} type="range" value={transform.order} min="0.25" max="3" step="0.25"/>
      <br/>
      <label htmlFor='transform_factor'>factor: {transform.factor} </label>
      <br/>
      <input id='transform_factor' onChange={(e)=>onChangeTransform({factor:e.target.value})} type="range" value={transform.factor} min="-1" max="1" step="0.1"/>
    </div>
  )
};
