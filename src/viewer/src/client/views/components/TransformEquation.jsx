import React, { Component } from 'react'
import { connect } from 'react-redux'
import styles from './Caption.scss';
import { getTransformEquation } from '../reducers/plotData'
import { format as d3Format } from 'd3-format'

class Equation extends Component {
  constructor(props) {
    super(props);
  }

  render(){
    let y = <i>y</i>
    if (this.props.yType == 'scaleLog'){
      y = <span> log(<i>y</i>)</span>
    }
    let x = <i>x</i>
    if (this.props.xType == 'scaleLog'){
      x = <span> log(<i>x</i>)</span>
    }
    return (
      <span>
        {y} = {y} {this.props.factor >= 0 ? '+' : '-'} {d3Format(",.4r")(Math.abs(this.props.factor))}{x} {this.props.intercept >= 0 ? '+' : '-'} {d3Format(",.4r")(Math.abs(this.props.intercept))}
      </span>
    )
  }
}

class TransformEquation extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return getTransformEquation(state)
    }
    this.mapDispatchToProps = dispatch => {
      return {

      }
    }
  }

  render(){
    const ConnectedEquation = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(Equation)
    return (
      <ConnectedEquation {...this.props}/>
    )
  }
}

export default TransformEquation
