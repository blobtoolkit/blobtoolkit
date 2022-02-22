import React from 'react';
import { connect } from 'react-redux'
import { getReferenceKites} from '../reducers/summary'
import { getShowReference} from '../reducers/reference'

export default class PlotRefKitesSVG extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => ({
        refs:getReferenceKites(state),
        showReference:getShowReference(state)
      })
    }
  }

  render(){
    const ConnectedRefKites = connect(
      this.mapStateToProps
    )(RefKitesSVG)
    return (
      <ConnectedRefKites {...this.props}/>
    )
  }
}

const RefKitesSVG = ({ refs, showReference }) => {
  if (!refs){
    return null
  }
  if (!showReference){
    return null
  }
  let refKites = []
  refs.coords.forEach((ref,i)=>{
    if (ref.x){
      let kite = <g key={i}
                    stroke={'#ccc'}
                    strokeWidth={1.5}
                    strokeDasharray={'3 6'}
                    strokeLinecap='round'
                    fill="#cccccc55">
                    <defs>
                      <path id={refs.ids[i]} d={'M'+ref.poly.map(c=>c[0]+' '+c[1]).join('L')} />
                    </defs>
                    <text fill='#999'
                          stroke='none'
                          fontSize={16}
                          textAnchor='start'
                          transform={'translate(4,-4)'}>
                      <textPath href={'#'+refs.ids[i]}>{refs.ids[i]}</textPath>
                    </text>
                    <line key={`${i}_x`}
                          x1={ref.x[0][0]}
                          y1={ref.x[0][1]}
                          x2={ref.x[1][0]}
                          y2={ref.x[1][1]}/>
                    <line key={`${i}_y`}
                          x1={ref.y[0][0]}
                          y1={ref.y[0][1]}
                          x2={ref.y[1][0]}
                          y2={ref.y[1][1]}/>
                    <polygon key={`${i}_poly`}
                        strokeWidth={1.5}
                        points={ref.poly.map(c=>c[0]+','+c[1]).join(' ')}/>
                  </g>
      refKites.push( kite )
    }
  })
  return (
    <g>
      {refKites}
    </g>
  )
}
