import React from 'react';
import { connect } from 'react-redux'
import styles from './Plot.scss'
import { format as d3format} from 'd3-format'
import Spinner from './Spinner'
import { plotPaths, plotText } from '../reducers/plotStyles'

export default class SnailPlotScale extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        {
          plotText: plotText(state),
          plotPaths: plotPaths(state)
        }
      )
    }
  }

  render(){
    const ConnectedSnailPlotScale = connect(
      this.mapStateToProps
    )(SnailPlotScaleComponent)
    return (
      <ConnectedSnailPlotScale {...this.props}/>
    )
  }
}

export const PlotTextInput = ({x,y,height,width,display,content,onChange,onKeyUp}) => {
  if (!display) return null
  return (
    <g>
      <foreignObject width={width} height={height}>
        <div xmlns="http://www.w3.org/1999/xhtml">
          <input
            className={styles.scale}
            type='number'
            value={content}
            onChange={e=>onChange(e.target.value)}
            onKeyUp={e=>onKeyUp(e.key,e.target.value)}
            style={{width,marginLeft:x,marginTop:y,padding:0,border:'3px solid'}}
          />
        </div>
      </foreignObject>
    </g>
  )

}
export class SnailPlotScaleComponent extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      circumference:props.scale.circumference,
      radius:props.scale.radius,
      editCircumference:false,
      editRadius:false
    }
  }

  render(){

    let items = []
    let title = this.props.title
    let scale = this.props.scale
    let offset = 15
    let w = 25
    let h = 25
    let gap = 10
    let ds = (
      <g transform={'translate('+0+','+0+')'}>
        <text style={this.props.plotText.snailLegendTitle}>{title}</text>
      </g>
    )
    let format = d3format(".2s")

    items.push(
      <g key={offset} transform={'translate(0,'+offset+')'}>
        <circle cx={w/2} cy={w/2} r={w/2} style={this.props.plotPaths.axis} fill='none' stroke='black'/>
        <line x1={w/2} y1={0} x2={w/2} y2={w/2} style={this.props.plotPaths.axis} fill='none' stroke='black'/>
        <rect
          x={w+gap}
          y={0}
          width='80px'
          height={w}
          style={Object.assign({},this.props.plotPaths.axis,{fill:'rgba(255,255,255,0)',stroke:'black',pointerEvents:'auto',cursor:'pointer'})}
          onClick={()=>this.setState({editCircumference:true})}
        />
      <text x={w+3*gap/2} y={w-gap/2} style={this.props.plotText.snailLegend}>{format(scale.circumference)}</text>
        <PlotTextInput
          x={w+gap}
          y={-1}
          height={w}
          width='120px'
          display={this.state.editCircumference}
          content={this.state.circumference}
          onChange={circumference=>this.setState({circumference})}
          onKeyUp={
            (key,value)=>{
              if (key === 'Enter') {
                this.props.onChangeCircumference(value)
                this.setState({editCircumference:false})
              }
            }
          }
        />
      </g>
    )
    offset += h + gap

    items.push(
      <g key={offset} transform={'translate(0,'+offset+')'}>
        <circle cx={w/2} cy={w/2} r={w/2} style={this.props.plotPaths.axis} fill='none' stroke='black'/>
        <line x1={w/2} y1={0} x2={w/2} y2={w/2} style={this.props.plotPaths.axis} fill='none' stroke='black'/>
        <rect
          x={w+gap}
          y={0}
          width='80px'
          height={w}
          style={Object.assign({},this.props.plotPaths.axis,{fill:'rgba(255,255,255,0)',stroke:'black',pointerEvents:'auto',cursor:'pointer'})}
          onClick={()=>this.setState({editRadius:true})}
        />
        <text x={w+3*gap/2} y={w-gap/2} style={this.props.plotText.snailLegend}>{format(scale.radius)}</text>
        <PlotTextInput
          x={w+gap}
          y={-1}
          height={w}
          width='120px'
          display={this.state.editRadius}
          content={this.state.radius}
          onChange={radius=>this.setState({radius})}
          onKeyUp={
            (key,value)=>{
              if (key === 'Enter') {
                this.props.onChangeRadius(value)
                this.setState({editRadius:false})
              }
            }
          }
        />
      </g>
    )
    offset += h + gap

    return (
      <g>
        {ds}
        {items}
      </g>
    )
  }
}
