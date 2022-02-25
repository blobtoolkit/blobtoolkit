import React from 'react';
import { connect } from 'react-redux'
import styles from './Plot.scss'

import { getCirclePlotDataForCategoryIndex }  from '../reducers/plotData'
import { getPlotResolution, getZScale }  from '../reducers/plotParameters'

export default class PlotBubblesCanvas extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        let val = getCirclePlotDataForCategoryIndex(state,this.props.x0)
        let obj = {}
        obj.color = val.color
        obj.data = val.data
        obj.scale = getZScale(state)
        obj.res = getPlotResolution(state)
        return obj
      }
    }
  }

  render(){
    const ConnectedBubblesCanvas = connect(
      this.mapStateToProps
    )(BubblesCanvas)
    return (
      <ConnectedBubblesCanvas {...this.props}/>
    )
  }
}

class BubblesCanvas extends React.Component {
  constructor(props) {
    super(props)
    let isChrome = /Chrome/.test(navigator.userAgent) && /Google Inc/.test(navigator.vendor)
    let isSafari = /Safari/.test(navigator.userAgent) && /Apple Computer/.test(navigator.vendor)
    if (isSafari){
      this.state = { webkit:true,height:900,width:900 }
    }
    else {
      this.state = { webkit:false,height:900,width:900 }
    }
    this.updateDimensions = this.updateDimensions.bind(this);
  }

  componentDidMount() {
    if (this.state.webkit){
      window.addEventListener("resize", this.updateDimensions);
    }
    this.updateDimensions()
  }
  componentDidUpdate() {
    this.updateCanvas()
  }

  updateDimensions(){
    if (this.state.webkit){
      let parent = this.canvas.parentNode.parentNode.parentNode.parentNode.parentNode.parentNode.parentNode.parentNode
      let width = parent.clientWidth * 900 / 1420
      let height = parent.clientHeight * 900 / 1420
      width = Math.min(width,height)
      height = width
      this.setState({height,width})
    }
    this.updateCanvas()
  }


  componentWillUnmount() {
    if (this.state.webkit){
      window.removeEventListener("resize", this.updateDimensions);
    }
  }

  updateCanvas() {
    let width = this.canvas.width
    let height = this.canvas.height
    const ctx = this.canvas.getContext('2d')
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.scale(1/window.devicePixelRatio,1/window.devicePixelRatio)
    ctx.clearRect(0, 0, width+100, height+100);
    ctx.globalAlpha=0.4
    ctx.fillStyle = this.props.color;
    ctx.lineWidth = 2;
    ctx.strokeStyle = 'rgb(89, 101, 111)';
    this.props.data.map(bubble => {
      ctx.beginPath();
      ctx.arc(bubble.x*width/900, bubble.y*width/900, bubble.r*width/900+1, 0, 2 * Math.PI, false);
      ctx.fill();
      ctx.stroke();
    })
  }


  render() {
    let scale = this.state.webkit ? this.state.width/(900) : this.state.width/900
    return (
      <canvas
        className={styles.main_canvas}
        ref={(elem) => { this.canvas = elem; }}
        width={this.canvas ? this.canvas.width : 900*window.devicePixelRatio }
        height={this.canvas ? this.canvas.height : 900*window.devicePixelRatio}
        style={
          {
            transform:'scale('+scale+')',
            transformOrigin: 'left top'
          }
        }
      />
    );
  }
}
