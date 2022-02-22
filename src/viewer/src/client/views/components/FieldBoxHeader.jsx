import React from 'react'
import styles from './Fields.scss';
import FieldBoxHeaderButton from './FieldBoxHeaderButton'
import FieldBoxHeaderFilterButtons from './FieldBoxHeaderFilterButtons'
import SVGIcon from './SVGIcon'
import xAxisIcon from './svg/xAxis.svg';
import yAxisIcon from './svg/yAxis.svg';
import zAxisIcon from './svg/zAxis.svg';
import categoryIcon from './svg/category.svg';
import sumActiveIcon from './svg/sumActive.svg';
import showIcon from './svg/show.svg';
import hideIcon from './svg/showHide.svg';
import selectInverseIcon from './svg/invert.svg';
import selectAllIcon from './svg/selectAll.svg';
import selectNoneIcon from './svg/selectNone.svg';
import invertSelectionIcon from './svg/invertSelection.svg';

class FieldBoxHeader extends React.Component {
  checkAxisStatus(axis){
    if (this.props.plot.axes &&
      this.props.plot.axes[axis] == this.props.fieldId){
        return true
      }
    return false
  }
  render(){
    let buttons
    let active = this.props.active && !this.props.isStatic
    if (active || this.props.type == 'selection'){
      let icons
      if (this.props.type == 'variable'){
        icons = (
          <span>
            <SVGIcon sprite={xAxisIcon} active={this.checkAxisStatus('x')} onIconClick={()=>{this.props.onAxisButtonClick('x',this.props.fieldId)}}/>
            <SVGIcon sprite={yAxisIcon} active={this.checkAxisStatus('y')} onIconClick={()=>{this.props.onAxisButtonClick('y',this.props.fieldId)}}/>
            <SVGIcon sprite={zAxisIcon} active={this.checkAxisStatus('z')} onIconClick={()=>{this.props.onAxisButtonClick('z',this.props.fieldId)}}/>
            { this.props.meta.parent != this.props.datasetId &&
              !this.props.clonedFrom &&
              <SVGIcon sprite={sumActiveIcon} onIconClick={()=>{this.props.onSumButtonClick(this.props.fieldId)}}/>
            }
          </span>
        )
      }
      else if (this.props.type == 'category'){
        icons = (
          <span>
            <SVGIcon sprite={categoryIcon} active={this.checkAxisStatus('cat')} onIconClick={()=>{this.props.onAxisButtonClick('cat',this.props.fieldId)}}/>
          </span>
        )
      }
      if (this.props.type == 'selection'){
        icons = (
          <span>
          {this.props.hideSelection ? '' :
            (<span>
              <SVGIcon sprite={invertSelectionIcon} onIconClick={()=>{this.props.invertSelection()}}/>
              <SVGIcon sprite={selectNoneIcon} onIconClick={()=>{this.props.selectNone()}}/>
              <SVGIcon sprite={selectAllIcon} onIconClick={()=>{this.props.selectAll()}}/>
            </span>)}
            <SVGIcon sprite={this.props.hideSelection ? hideIcon : showIcon} active={this.props.hideSelection} onIconClick={()=>{this.props.toggleSelection(!this.props.hideSelection)}}/>
          </span>
        )
      }
      buttons = (
        <div className={active ? styles.header_buttons : styles.alt_header_buttons}>
          {icons}
          <FieldBoxHeaderFilterButtons {...this.props}/>
        </div>
      )
    }
    else {
      buttons = (
        <div className={styles.header_type}
          data-tip data-for='header-type'
          >{this.props.type}</div>
        )
    }
    return (
      <div className={styles.header + ' ' + (active ? styles.active : '')}>
        <span className={styles.left_half}
          data-tip data-for={(this.props.active ? 'active-' : '') + 'field-header'}
          onClick={()=>{this.props.onHeaderClick()}}>
          <h1 className={styles.inline}>{this.props.fieldId}</h1>
        </span>
        {buttons}
      </div>
    );
  }
}

export default FieldBoxHeader
