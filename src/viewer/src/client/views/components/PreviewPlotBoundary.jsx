import React from 'react';
import { connect } from 'react-redux'
import { getMainPlotData }  from '../reducers/plotData'
import { plotText }  from '../reducers/plotStyles'
import styles from './Plot.scss'
import {Axis, axisPropsFromTickScale, LEFT, RIGHT, TOP, BOTTOM} from 'react-d3-axis';
import { getPreviewDataForFieldId } from '../reducers/field'
import { format as d3format} from 'd3-format'

export default class PreviewPlotBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      let preview = getPreviewDataForFieldId(state,props.xLabel) || {}
      return (state, props) => (
        {
          ...preview,
          plotText: plotText(state)
        }
      )
    }
  }

  render(){
    const PlotBoundary = connect(
      this.mapStateToProps
    )(PlotOutline)
    return (
      <PlotBoundary {...this.props}/>
    )
  }
}

class PlotOutline extends React.Component {
  render(){
    if (!this.props.bars) return (<g></g>)
    let yScale = this.props.yScale.copy()
    let height = this.props.dimensions.height
    let width = this.props.dimensions.width
    yScale.range([height,0])
    yScale.domain([0,this.props.max])
    let fontSize = this.props.plotText.axisTick.fontSize.replace('px','')/2.5
    let format = d3format(".2s")
    let altFormat = d3format(".2f")
    let xScale = this.props.details.xScale.copy()
    let ticks, labels
    let clamp
    if (this.props.details.type == 'variable'){
      let count = 25
      xScale.range([0,count])
      let thresh = Array.from(Array(count-1).keys()).map((n)=>{let v = xScale.invert((n+1)); return v > 0.001 && v < 1 ? altFormat(v) : format(v)});
      ticks = this.props.bars.map((bar,i)=>{
        let x = bar.x
        let line = <line key={i} stroke='black' x1={x} x2={x} y2={6}/>
        if (this.props.details.meta.clamp &&
          this.props.details.meta.clamp > this.props.details.meta.limit[0]
          && i < 2) line = null
        return line
      })
      let w = width / count
      labels = thresh.map((t,i)=>{
        let text = <text transform={'translate('+(w+i*w)+',10),rotate(90)'} key={i} textAnchor='start' dominantBaseline='middle' style={{fontSize}}>{t}</text>
        if (this.props.details.meta.clamp &&
          this.props.details.meta.clamp > this.props.details.meta.limit[0]
          && i < 1) text = null
        return text
      })
      ticks = ticks.slice(1)
      clamp = this.props.details.meta.clamp &&
        this.props.details.meta.clamp > this.props.details.meta.limit[0] && (
        <g stroke='black' transform='translate(0,3)'>
          <line x1={this.props.bars[0].x} x2={this.props.bars[0].x} y2={6}/>
          <line x1={this.props.bars[0].x} x2={this.props.bars[1].x}/>
          <line x1={this.props.bars[1].x} x2={this.props.bars[1].x} y2={6}/>
          <text stroke='none' transform={'translate('+(this.props.bars[1].x/2)+',10)'} textAnchor='middle' dominantBaseline='hanging' style={{fontSize}}>&lt; {this.props.details.meta.clamp}</text>
        </g>
      )
    }
    else if (this.props.details.type == 'category'){
      let count = 10
      xScale.range([0,count])
      let thresh = Array.from(Array(count-1).keys()).map((n)=>{let v = xScale.invert((n+1)); return v > 0.001 && v < 1 ? altFormat(v) : format(v)});
      ticks = this.props.bars.map((bar,i)=>{
        let x = bar.x
        return (
          <line key={i} stroke='black' x1={x} x2={x} y2={6}/>
        )
      })
      let w = width / count
      labels = this.props.bins.map((bin,i)=>{
        return (
          <text transform={'translate('+(w/2+i*w)+',5),rotate(90)'} key={i} textAnchor='start' dominantBaseline='middle' style={{fontSize}}>{bin.id}</text>
        )
      })
      ticks = ticks.slice(1)
    }

    return (
      <g>
        <g  transform={'translate(-5)'} >
          <Axis {...axisPropsFromTickScale(yScale, 4)} style={{orient: LEFT, tickFontSize: fontSize}} format={format}/>
        </g>
        <g  transform={'translate('+(width+5)+')'} >
          <Axis {...axisPropsFromTickScale(yScale, 4)} style={{orient: RIGHT, tickFontSize: 0}}/>
        </g>
        <g  transform={'translate(0,'+(height+5)+')'} >
          <line stroke='black' x2={width}/>
          {ticks}
          {labels}
          {clamp}
          <text transform={'translate('+(width/2)+',65)'} textAnchor='middle'>{this.props.xLabel}</text>
        </g>

      </g>
    )
  }
}
