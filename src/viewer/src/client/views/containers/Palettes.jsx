import React from 'react'
import { connect } from 'react-redux'
import { getAllPalettes, selectPalette, editPalette, choosePalette, chooseColors, getSelectedPalette, getUserPalette } from '../reducers/color'
import PalettesComp from '../components/Palettes'

class Palettes extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        let palettes = getAllPalettes(state)
        let selected = getSelectedPalette(state)
        let userPalette = getUserPalette(state)
        return {...palettes,selected,userPalette}
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        selectPalette: (id) => {
          return dispatch(choosePalette(id))
        },
        editPalette: (obj) => {return dispatch(chooseColors(obj))}
      }
    }
  }

  render(){
    const ConnectedPalettes = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(PalettesComp)
    return (
      <ConnectedPalettes name={'palettes'} {...this.props}/>
    )
  }
}

export default Palettes
