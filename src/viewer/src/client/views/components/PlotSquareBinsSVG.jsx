import React from 'react';
import { connect } from 'react-redux'
import PlotSquareBinSVG from './PlotSquareBinSVG'
import { getScatterPlotDataBySquareBinByCategory }  from '../reducers/plotSquareBins'

export default class PlotSquareBinsSVG extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        getScatterPlotDataBySquareBinByCategory(state)
      )
    }
  }

  render(){
    const ConnectedSquareBins = connect(
      this.mapStateToProps
    )(SquareBinsSVG)
    return (
      <ConnectedSquareBins {...this.props}/>
    )
  }
}

const SquareBinsSVG = ({ data = [], css = '' }) => {
  return (
    <g transform='translate(50, 50)'>
    {data.map(square =>
      <PlotSquareBinSVG key={square.id} {...square} />
    )}
    </g>
  )
}
