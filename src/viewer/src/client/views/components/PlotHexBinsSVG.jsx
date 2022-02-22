import React from 'react';
import { connect } from 'react-redux'
import PlotHexBinSVG from './PlotHexBinSVG'
import { getScatterPlotDataByHexBinByCategory }  from '../reducers/plotHexBins'

export default class PlotHexBinsSVG extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        //console.log(getScatterPlotDataByHexBin(state))
        getScatterPlotDataByHexBinByCategory(state)
        // plotGraphics:getPlotGraphics(state)
      )
    }
  }

  render(){
    const ConnectedHexBins = connect(
      this.mapStateToProps
    )(HexBinsSVG)
    return (
      <ConnectedHexBins {...this.props}/>
    )
  }
}

const HexBinsSVG = ({ data = [], css = '' }) => (
    <g transform='translate(50, 50)'>
      {data.map(hex =>
        <PlotHexBinSVG key={hex.id} {...hex} css={css} />
      )}
    </g>
);
