import React from 'react'
import styles from './Palette.scss'
import { SketchPicker } from 'react-color'
import { queryToStore } from '../querySync'


class PaletteSwatch extends React.Component {
  constructor(props) {
    super(props);
    let color = this.props.color.replace(/[rgba\(\)]/g,'').split(',');
    this.state = {
      displayColorPicker: false,
      color: {
        r: color[0],
        g: color[1],
        b: color[2],
        a: color[3] || 1
      },
    };
  }

  componentDidUpdate(nextProps) {
    if (this.props.color != nextProps.color && !this.state.displayColorPicker){
      let color = this.props.color.replace(/[rgba\(\)]/g,'').split(',');
      this.setState({
        displayColorPicker: false,
        color: {
          r: color[0],
          g: color[1],
          b: color[2],
          a: color[3] || 1
        },
      })
    }

  }

  handleClick() {
    this.setState({ displayColorPicker: !this.state.displayColorPicker })
  };

  handleClose() {
    let rgb = 'rgb('+this.state.color.r+','+this.state.color.g+','+this.state.color.b+')'
    this.setState({ displayColorPicker: false })
  };

  render(){
    const presetColors = this.props.defaultColors.concat([
      "rgb(255,255,0)",
      "rgb(0,0,0)",
      "rgb(63,63,63)",
      "rgb(127,127,127)",
      "rgb(191,191,191)",
      "rgb(255,255,255)"
    ])
    const handleChange = (color) => {
      this.setState({ color: color.rgb })
    };
    const handleChangeComplete = (color) => {
      let rgba = 'rgba('+color.rgb.r+','+color.rgb.g+','+color.rgb.b+','+color.rgb.a+')'
      this.props.editPalette(rgba)
    };
    return (
      <div className={styles.block}
        style={{width:this.props.width}}
        data-tip data-for='edit-swatch'>
        <div className={ styles.swatch } onClick={ ()=>{this.handleClick()} }>
          <div className={ styles.color }
            style={{backgroundColor:`rgba(${ this.state.color.r }, ${ this.state.color.g }, ${ this.state.color.b }, ${ this.state.color.a })`}} />
        </div>
        { this.state.displayColorPicker ? <div className={ styles.popover }>
          <div className={ styles.cover } onClick={ ()=>{this.handleClose()} }/>
          <SketchPicker color={ this.state.color } presetColors={presetColors} onChange={ handleChange } onChangeComplete={ handleChangeComplete } />
        </div> : null }

      </div>
    )
  }
}

export default PaletteSwatch
