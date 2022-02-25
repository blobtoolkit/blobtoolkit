import React from 'react';
import { connect } from 'react-redux'
import styles from './Layout.scss'

class NumericInput extends React.Component {
  constructor(props) {
    super(props);
    let value = isNaN(this.props.initialValue) ? 0 : this.props.initialValue
    this.state = {value}
  }

  componentDidUpdate(nextProps) {
    let value = isNaN(this.props.initialValue) ? 0 : this.props.initialValue
    let npvalue = isNaN(nextProps.initialValue) ? 0 : nextProps.initialValue
    if (npvalue != value){
      this.setState({value})
    }
  }
  render(){
    let minValue = this.props.minValue || 0
    return (
      <input
        type='number'
        className={styles.range_input}
        value={this.state.value}
        onChange={e=>{
          this.setState({value:e.target.value})
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
            this.props.onChange(e.target.value, minValue)
          }
        }
      />
    )
  }
}

export default NumericInput
