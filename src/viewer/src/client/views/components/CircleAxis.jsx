import React from 'react';
import { connect } from 'react-redux'
import { scaleLinear as d3ScaleLinear } from 'd3-scale';
import { radialLine as d3RadialLine } from 'd3-shape';
import { arc as d3Arc } from 'd3-shape';
import { format as d3Format } from 'd3-format'
import { plotPaths, plotText } from '../reducers/plotStyles'


export class CircleAxis extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => (
        {
          plotText: plotText(state),
          plotPaths: plotPaths(state)
        }
      )
    }
  }

  render(){
    const ConnectedCircleAxis = connect(
      this.mapStateToProps
    )(CircleAxisComponent)
    return (
      <ConnectedCircleAxis {...this.props}/>
    )
  }
}

export const CircleAxisComponent = ({domain=[0, 100], range=[0, 360], unit='%', major=10, minor=2, radius=100, inner=false, outer=true, plotPaths, plotText}) => {
  let cScale = d3ScaleLinear().domain(domain).range(range.map(x=>x*2*Math.PI/360))
  let rScale = d3ScaleLinear().domain([0, 100]).range([0, radius])
  let circ = d3Arc()({
    startAngle: cScale(range[0]),
    endAngle: cScale(range[1]),
    innerRadius: rScale(0),
    outerRadius: rScale(100)
  })
  let start = d3RadialLine()(
    [
      [cScale(domain[0]),rScale(0)],
      [cScale(domain[0]),rScale(100)]
    ]
  )
  let ticks = []
  let labels = []
  let rmin = inner ? 94 : 100
  let rmax = outer ? 106 : 100
  for (let i = domain[0]; i <= domain[1]; i += major){
    ticks.push(d3RadialLine()(
      [
        [cScale(i),rScale(rmin)],
        [cScale(i),rScale(rmax)]
      ]
    ))
    let text = i
    if (i == domain[0]){
      text += unit
    }
    if (i == domain[1]){
      text = ''
    }
    labels.push(
      {
        path: d3RadialLine()(
          ([...Array(10).keys()]).map(j=>(
            [cScale(i+j-9),rScale(rmax+3)],
            [cScale(i+j-9),rScale(rmax+3)]
          ))
        ),
        text: text,
        fontSize: rScale(plotText.axisTick.fontSize.replace('px',''))
      }
    )
  }
  rmin = inner ? 97 : 100
  rmax = outer ? 103 : 100
  let minorTicks = []
  for (let i = domain[0]; i <= domain[1]; i += minor){
    if (i % major != 0){
      minorTicks.push(d3RadialLine()(
        [
          [cScale(i),rScale(rmin)],
          [cScale(i),rScale(rmax)]
        ]
      ))
    }
  }
  return (
    <g>
      <path d={circ} style={plotPaths.axis} stroke='black'/>
      <path d={start} style={plotPaths.axis} stroke='black'/>
      {
        ticks.map((d,i)=>(
          <path d={d} key={i} style={plotPaths.axis} stroke='black'/>
        ))
      }
      {
        minorTicks.map((d,i)=>(
          <path d={d} key={i} style={plotPaths.fine} stroke='black'/>
        ))
      }
      {
        labels.map((d,i)=>{
          return (
            <g key={i}>
              <path d={d.path} id={'path_'+i} style={plotPaths.axis} stroke='none'/>
              <text style={Object.assign({}, plotText.axisLabel, {fontSize:d.fontSize})}>
                <textPath
                      xlinkHref={'#path_'+i}
                      textAnchor={'end'}
                      startOffset={'100%'}>
                  {d.text}
                </textPath>
              </text>
            </g>
          )
        })

      }

    </g>
  )
}
