import React from 'react'
import ReactDOM from 'react-dom'
import styles from './Filters.scss';
// import { DraggableCore } from 'react-draggable'
import Pointable from 'react-pointable'
import { format as d3format } from "d3-format";

class FilterDisplayRange extends React.Component {
  render(){
    return (
      <FilterHandles {...this.props}>
      </FilterHandles>
    )
  }
}

class FilterHandles extends React.Component {
  render() {
    return (
      <div className={styles.handles_container} ref='handleDiv'>
        <FilterHandle {...this.props} key='right' handlePosition='right' />
        <FilterHandle {...this.props} key='left' handlePosition='left' />
      </div>
    )
  }
}

class FilterHandle extends React.Component {
  constructor(props) {
    super(props);
    this.state = {offsetX:0,mouseDown:false,css:''}
  }
  setMouseDown(bool){
    this.setState({mouseDown:bool,css:bool?styles.wide_handle:''})
  }
  bound(){
    return this.props.handlePosition == 'right' ? 1 : 0
  }
  boundPx(){
    return this.props.xScale(this.props.filterRange[this.bound()]);
  }
  updateRange(value,bound) {
    let range = this.props.filterRange.slice(0);
    this.setState({offsetX:0})
    let index = bound + 1
    let v = value
    if (bound == 0){
      v = Math.max(this.props.filterLimit[0],value)

    }
    else {
      v = Math.min(this.props.filterLimit[1],value)
    }

    if (v != value){
      index *= -1
    }
    if (bound == 0 && this.props.meta.clamp){
      let clamp = this.props.meta.clamp // || this.props.filterLimit[0]
      if (value < clamp){
        if (value > range[0]){
          v = clamp
        }
        else {
          v = this.props.filterLimit[0]
          index = -1
        }
      }
    }
    range[bound] = d3format(".3r")(v)
    this.props.onUpdateRange(this.props.filterId,range,index)
  }
  render(){
    return (
      <Pointable
        tagName='span'
        onPointerMove={(e)=>{
          e.preventDefault()
          if (this.state.mouseDown){
            this.setState({offsetX:e.x-this.boundPx()})
          }
        }}
        onPointerDown={(e)=>{
          e.preventDefault()
          this.setMouseDown(true)
        }}
        onPointerUp={(e)=>{
          e.preventDefault()
          if (this.state.mouseDown){
            this.setMouseDown(false)
            this.updateRange(this.props.xScale.invert(e.x),this.bound())
          }
        }}
        onPointerLeave={(e)=>{
          e.preventDefault()
          if (this.state.mouseDown){
            this.setMouseDown(false)
            this.updateRange(this.props.xScale.invert(e.x),this.bound())
          }
        }}
        >
        <div style={{left: (this.boundPx()+this.state.offsetX)+'px'}}
        className={styles.handle+' '+styles[this.props.handlePosition]+' '+this.state.css}>
          <div className={styles.arrows} data-tip data-for='draggable-arrow'>
            &lt;&nbsp;&gt;
          </div>
        </div>
      </Pointable>
    )
  }
}


export default FilterDisplayRange
