import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'
import paletteStyles from './Palette.scss'
import MenuDisplaySimple from './MenuDisplaySimple'
import TextIcon from './TextIcon'
import { getTransformFunctionParams,
  chooseTransformFunction } from '../reducers/plotParameters'
import { getKitePlotData, scaleFactor } from '../reducers/plotData'
import { getColorPalette } from '../reducers/color'
import SVGIcon from './SVGIcon'
import kiteIcon from './svg/kiteSolid.svg'
import NumericInput from './NumericInput'
import { format as d3Format } from 'd3-format'

const dataset_table = DATASET_TABLE || false

class KiteMenu extends React.Component {
  constructor(props) {
    super(props);
  }

  updateTransform(factor, intercept){
    let origin = {x:0, y:this.props.kite.yScale(this.props.kite.yScale.domain()[1]-intercept)}
    if (intercept != 0){
      if (this.props.kite.yScale.type == 'scaleLog'){
        intercept = this.props.kite.yScale(10**intercept)
      }
      else {
        intercept = this.props.kite.yScale(intercept)
      }
    }
    factor = factor == 0 ? 0 : factor /= this.props.scaleFactor.factor
    intercept = d3Format(",.4r")(intercept)
    factor = d3Format(",.4r")(factor)
    this.props.onSelectKite({factor, intercept, origin})
  }


  render(){
    let {params, kite, palette, scaleFactor} = this.props
    if (!kite || !kite.equations) return null
    let buttons = []
    let colors = palette.colors
    kite.equations.forEach((equation,i)=>{
      if (equation.hasOwnProperty('intercept')){
        let eqn = Object.assign({},equation)
        if (Math.abs(eqn.factor) > 0.0001){
          eqn.factor += params.factor
          eqn.intercept += params.intercept
        }
        else {
          eqn.factor = params.factor
          eqn.intercept = params.intercept
        }
        buttons.push(
          <span key={i}
                style={{color:colors[i]}}>
            <SVGIcon sprite={kiteIcon}
                     active={i == params.index}
                     invert={true}
                     onIconClick={()=>this.props.onSelectKite(eqn)}/>
                 </span>

        )
      }
    })
    let reset = params.factor != 0 || params.intercept != 0
    let intercept
    if (params.intercept == 0){
      intercept = 0
    }
    else if (kite.yScale.type == 'scaleLog') {
      intercept = Math.log10(kite.yScale.invert(params.intercept))
    }
    else {
      intercept = kite.yScale.invert(params.intercept)
    }
    kite.yScale.invert(params.intercept)
    let factor
    if (params.factor == 0){
      factor = 0
    }
    else {
      factor = params.factor * scaleFactor.factor
    }
    return (
      <div style={{marginBottom:'0.5em'}}>
        <MenuDisplaySimple name='transform' key={1} invert={false}>
          <div className={styles.full_height}>
            {reset && <span style={{position:'absolute', top:'2em', left:'1em', zIndex:5}}>
                        <TextIcon invert={true}
                                  title={'reset'}
                                  active={false}
                                  onIconClick={()=>this.props.onSelectKite({})}/>
                              </span>}
            {kite.yScale.type == 'scaleLog' && <span>log(<i>y</i>)</span> || <i>y</i>} = <NumericInput initialValue={d3Format(",.4r")(factor)} onChange={(value)=>this.updateTransform(value, intercept)}/>
            {kite.xScale.type == 'scaleLog' && <span>log(<i>x</i>)</span> || <i>x</i>} + <NumericInput initialValue={d3Format(",.4r")(intercept)} onChange={(value)=>this.updateTransform(factor, value)}/>
          </div>
        </MenuDisplaySimple>
        <MenuDisplaySimple key={2} invert={true}>
          {buttons}
        </MenuDisplaySimple>
      </div>
    )
  }
}

class MenuDisplayKite extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        params: getTransformFunctionParams(state),
        kite: getKitePlotData(state),
        palette: getColorPalette(state),
        scaleFactor: scaleFactor(state)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        onSelectKite: obj => {
          if (obj.hasOwnProperty('intercept')) obj.intercept = d3Format(",.3f")(obj.intercept)
          if (obj.hasOwnProperty('factor')) obj.factor = d3Format(",.4r")(obj.factor)
          dispatch(chooseTransformFunction(obj))
        }
      }
    }
  }

  render(){
    const ConnectedKiteMenu = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(KiteMenu)
    return <ConnectedKiteMenu {...this.props}/>
  }
}

export default MenuDisplayKite
