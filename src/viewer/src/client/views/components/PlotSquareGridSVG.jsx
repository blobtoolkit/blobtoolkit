import React from 'react';
import { connect } from 'react-redux'
import styles from './Plot.scss'
import { getSelectedSquareGrid } from '../reducers/plotSquareBins'
import { addRecords, removeRecords, setSelectSource } from '../reducers/select'
import {
  grid,
  gridShape,
  gridShapeSelected,
  gridShapePartSelected
} from '../reducers/plotStyles'

export default class PlotSquareGridSVG extends React.Component {
  constructor(props) {
    super(props);
    this.state = {mouseDown:false,addRecords:true}
    this.mapStateToProps = () => {
      return (state, props) => (
        {
          squareGrid: getSelectedSquareGrid(state),
          gridStyle: grid(state),
          gridShape: gridShape(state),
          gridShapeSelected: gridShapeSelected(state),
          gridShapePartSelected: gridShapePartSelected(state)
        }
      )
    }
    this.mapDispatchToProps = dispatch => {
      return {
        onClickCell:(arr) => {
          if (this.state.addRecords){
            dispatch(setSelectSource('square'))
            dispatch(addRecords(arr))
          }
          else {
            dispatch(setSelectSource('square'))
            dispatch(removeRecords(arr))
          }
        }
      }
    }
  }

  setMouseDown(bool){
    this.setState({mouseDown:bool})
  }

  setAddRecords(bool){
    this.setState({addRecords:bool})
  }

  render(){
    const ConnectedHexGrid = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(SquareGridSVG)
    return (
      <ConnectedHexGrid {...this.props}
        mouseDown={this.state.mouseDown} setMouseDown={(bool)=>this.setMouseDown(bool)}
        addRecords={this.state.addRecords} setAddRecords={(bool)=>this.setAddRecords(bool)}
        />
    )
  }
}
const SquareGridSVG = ({ squareGrid,
                      onClickCell,
                      mouseDown,
                      setMouseDown,
                      setAddRecords,
                      gridStyle,
                      gridShape,
                      gridShapeSelected,
                      gridShapePartSelected }) => {
  let squares = []
  let data = squareGrid.data
  data.forEach((datum,i)=>{
    let style = gridShape
    if (datum.selected){
      style = gridShapeSelected
      if (datum.selected < datum.ids.length){
        style = gridShapePartSelected
      }
    }
    squares.push(
      <rect key={i}
        style={style}
        x={datum.x}
        y={datum.y}
        height={datum.height}
        width={datum.width}
        onMouseOver={()=>{if (mouseDown){onClickCell(datum.ids)}}}
        onMouseDown={()=>{
          if (datum.selected){
            setAddRecords(false)
          }
          else {
            setAddRecords(true)
          }
          setMouseDown(true)
        }}
        onMouseUp={()=>{setMouseDown(false)}}
      />)
  })
  return (
    <g transform='translate(50, 50)' style={{gridStyle}}>
    {squares}
    </g>
  )
};
