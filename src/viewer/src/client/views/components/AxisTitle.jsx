import React from 'react'
import { connect } from 'react-redux'
import { getAxisTitle } from '../reducers/plotData'
import styles from './Plot.scss'
import { plotText } from '../reducers/plotStyles'

export default class AxisTitle extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        {
          plotText: plotText(state)
        }
      )
    }
  }

  render(){
    const ConnectedAxisTitle = connect(
      this.mapStateToProps
    )(AxisTitleComponent)
    return (
      <ConnectedAxisTitle {...this.props}/>
    )
  }
}

const AxisTitleComponent = ({ axis, title, side=1000, plotText }) => {
  let params = {}
  params.transform = axis == 'x' ? 'translate(500,1085)' : 'translate(-90,500),rotate(90)'
  return (
    <g {...params}>
      <text style={plotText.axisTitle} transform={'scale('+side/1000+')'}>
        {title}
      </text>
    </g>
  )
};
