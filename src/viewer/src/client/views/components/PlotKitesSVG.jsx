import React from 'react';
import { connect } from 'react-redux'
import { getKitePlotData }  from '../reducers/plotData'

export default class PlotKitesSVG extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        getKitePlotData(state)
      )
    }
  }

  render(){
    const ConnectedKites = connect(
      this.mapStateToProps
    )(KitesSVG)
    return (
      <ConnectedKites {...this.props}/>
    )
  }
}

const KitesSVG = ({ coords, bins, colors }) => {
  if (!coords){
    return null
  }
  let straight = []
  let paths = []
  coords.forEach((group,i)=>{
    if (coords[i].x){
      let kite = <g key={i}
                    style={{strokeWidth:"1px"}}
                    transform={`rotate(${coords[i].angle},${coords[i].y[0][0]},${coords[i].x[0][1]})`}
                    stroke={colors[i]}
                    fill="none">
                    <line key={`${i}_x`}
                          style={{strokeWidth:"3px"}}
                          x1={coords[i].x[0][0]}
                          y1={coords[i].x[0][1]}
                          x2={coords[i].x[1][0]}
                          y2={coords[i].x[1][1]}/>
                    <line key={`${i}_y`}
                          x1={coords[i].y[0][0]}
                          y1={coords[i].y[0][1]}
                          x2={coords[i].y[1][0]}
                          y2={coords[i].y[1][1]}/>
                    <polygon key={`${i}_poly`}
                             style={{strokeWidth:`3px`}}
                             points={coords[i].poly.map(c=>c[0]+','+c[1]).join(' ')}/>
                  </g>
      paths.push( kite )
      // straight.push( <g key={i}
      //                style={{strokeWidth:"1px",opacity:'0.25'}}
      //                >
      //               {kite}
      //             </g>)
    }
  })
  return (
    <g>
      {straight}
      {paths}
    </g>
  )
}
