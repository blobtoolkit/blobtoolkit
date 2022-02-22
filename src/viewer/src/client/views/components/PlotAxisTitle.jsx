import React from 'react'
import { connect } from 'react-redux'
import { getAxisTitle } from '../reducers/plotData'
import AxisTitle from './AxisTitle'

export default class PlotAxisTitle extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => ({
        title:getAxisTitle(state,props.axis)
      })
    }
  }

  render(){
    const ConnectedAxisTitle = connect(
      this.mapStateToProps
    )(AxisTitle)
    return (
      <ConnectedAxisTitle {...this.props} />
    )
  }
}
