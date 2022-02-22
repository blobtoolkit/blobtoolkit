import React from 'react'
import styles from './Filters.scss'
import { format as d3format } from "d3-format";

class FilterControlRange extends React.Component {
  constructor(props) {
    super(props);
    this.state = {low:props.filterRange[0],high:props.filterRange[1],clamp:props.meta.clamp}
  }
  componentDidUpdate(nextProps) {
    if (nextProps.filterRange[0] != this.props.filterRange[0] || nextProps.filterRange[1] != this.props.filterRange[1]){
      this.setState({low:nextProps.filterRange[0],high:nextProps.filterRange[1]})
    }
  }
  isNumeric(n) {
    if ((typeof n === 'undefined') || n == 'NaN') return false
    return !isNaN(parseFloat(n)) && isFinite(n)
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

  resetInvert(){
    let values = {}
    values[this.props.filterId+'--Inv'] = false
    let remove = [
      this.props.filterId+'--Inv'
    ]
    remove.forEach(key=>{
      delete this.props.parsed[key]
    })
    this.props.onChangeAxisRange(values, remove)
    this.props.changeQueryParams(this.props.parsed)
  }

  resetLimits() {
    let values = {}
    values[this.props.filterId+'--LimitMin'] = this.props.range[0]
    values[this.props.filterId+'--Min'] = this.props.range[0]
    values[this.props.filterId+'--LimitMax'] = this.props.range[1]
    values[this.props.filterId+'--Max'] = this.props.range[1]
    let remove = [
      this.props.filterId+'--LimitMin',
      this.props.filterId+'--Min',
      this.props.filterId+'--LimitMax',
      this.props.filterId+'--Max'
    ]
    remove.forEach(key=>{
      delete this.props.parsed[key]
    })
    this.props.onChangeAxisRange(values, remove)
    this.props.changeQueryParams(this.props.parsed)
  }

  render() {
    let reset
    let resetBounds
    if (this.props.filterRange[0] != this.props.fieldLimit[0]
        || this.props.filterRange[1] != this.props.fieldLimit[1]
        || this.props.invert){
      reset = (<span className={styles.reset}
                onClick={()=>{
                  this.updateRange(this.props.fieldLimit[0]*1-1,0);
                  this.updateRange(this.props.fieldLimit[1]*1+1,1);
                  this.resetInvert();
                }}>
                 reset
               </span>)
    }
    else if (this.props.range[0] != this.props.fieldLimit[0] || this.props.range[1] != this.props.fieldLimit[1]){
      resetBounds = (<span className={styles.reset}
                onClick={()=>{this.resetLimits(this.props)}}>
                 reset range
               </span>)
    }
    let clamp
    if (this.props.meta.clamp && this.props.meta.clamp > 0 && this.props.meta.limit[0] < this.props.meta.clamp){
       clamp = (
         <input
           type='number'
           className={styles.clamp_input}
           value={this.state.clamp}
           onChange={e=>{
             let value = e.target.value;
             let clamp
             this.setState({clamp:value})
           }}
           onKeyUp={
             (e)=>{
               if (e.key === 'Enter') {
                 e.target.blur()
               }
             }
           }
           onBlur={
             (e)=>{
               let value = e.target.value;
               let clamp
               if (this.props.filterLimit[0] < value && this.props.filterLimit[1] > value){
                 this.setState({clamp:value})
                 this.props.onUpdateClamp(this.props.meta.id,value)
               }
               else {
                 this.props.onUpdateClamp(this.props.meta.id,this.props.meta.clamp)
               }
             }
           }
           data-tip data-for='clamp-input' />
       )
    }
    return (
      <div className={styles.headliner}>
        {clamp}
        <div rel={this.props.filterRange[0]} className={styles.range}>
          <input
            type='number'
            className={styles.range_input}
            value={this.state.low}
            onChange={e=>{
              this.setState({low:e.target.value})
            }}
            onKeyUp={
              (e)=>{
                if (e.key === 'Enter') {
                  e.target.blur()
                }
              }
            }
            onBlur={
              (e)=>{
                // let range = this.props.filterRange;
                let value = e.target.value;
                // let index = 1
                // if (this.props.filterLimit[0] <= value){
                //   range[0] = value
                // }
                // else {
                //   range[0] = this.props.filterLimit[0]
                //   this.setState({low:range[0]})
                //   index = -1
                // }
                this.updateRange(value,0)
              }
            }
            data-tip data-for='range-input' />
          &nbsp;:&nbsp;
          <input
            type='number'
            className={styles.range_input}
            value={this.state.high}
            onChange={e=>{
              this.setState({high:e.target.value})
            }}
            onKeyUp={
              (e)=>{
                if (e.key === 'Enter') {
                  e.target.blur()
                }
              }
            }
            onBlur={
              (e)=>{
                // let range = this.props.filterRange;
                let value = e.target.value;
                // console.log(value)
                // let index = 2
                // if (this.props.filterLimit[1] >= value){
                //   range[1] = value
                // }
                // else {
                //   range[1] = this.props.filterLimit[1]
                //   this.setState({high:range[1]})
                //   index = -2
                // }
                // this.props.onUpdateRange(this.props.filterId,range,index)
                this.updateRange(value,1)
              }
            }
            data-tip data-for='range-input' />
          {reset}
          {resetBounds}
        </div>
      </div>
    )
  }
}


export default FilterControlRange
