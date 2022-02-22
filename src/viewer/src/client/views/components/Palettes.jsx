import React from 'react'
import styles from './Menu.scss'
import Palette from './Palette'
import TextIcon from './TextIcon'

class Palettes extends React.Component {
  render(){
    let palettes = []
    this.props.allIds.forEach(id => {
      if (this.props.selected==id && id != 'default'){
        palettes.push(
          <Palette
            id={id}
            key={id}
            active={this.props.selected==id}
            location={this.props.location}
            colors={this.props.byId[id]}
            selectPalette={this.props.selectPalette}
            defaultColors={this.props.byId['default']}
            editPalette={(obj)=>this.props.editPalette(obj)}
            />
          )
      }
    })
    return (
      <div className={styles.outer}>
        <div className={styles.header}>
          <h1 className={styles.inline}>{this.props.name}</h1>
          <div style={{float:'right'}}>
            <TextIcon title='default' active={this.props.selected == 'default'} onIconClick={()=>this.props.selectPalette('default')}/>
            <TextIcon title='custom' active={this.props.selected == 'user'} onIconClick={()=>this.props.selectPalette('user')}/>
          </div>
        </div>
        <div className={styles.content}>
          {palettes}
        </div>
      </div>
    )
  }
}

export default Palettes
