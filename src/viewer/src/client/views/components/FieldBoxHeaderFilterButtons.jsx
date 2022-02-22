import React from 'react'
import { connect } from 'react-redux'
import styles from './Fields.scss';
import SVGIcon from './SVGIcon'
import invertIcon from './svg/invert.svg';
import { editFilter } from '../reducers/filter'
import { getDetailsForFilterId } from '../reducers/preview'
import { queryToStore } from '../querySync'

class FieldBoxHeaderFilterButtons extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        return getDetailsForFilterId(state, props.fieldId)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        toggleInvert: (id,value) => {
          value = !value
          if (value){
            dispatch(queryToStore({values:{[id+'--Inv']:value},action:'FILTER'}))
          }
          else {
            dispatch(queryToStore({values:{},remove:[id+'--Inv'],action:'FILTER'}))
          }
          dispatch(editFilter({id,invert:value}))
        }
      }
    }
  }

  render(){
    const ConnectedButtons = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(FilterButtons)
    return (
      <ConnectedButtons filterId={this.props.fieldId} {...this.props}/>
    )
  }
}

class FilterButtons extends React.Component {
  checkAxisStatus(axis){
    if (this.props.plot.axes &&
      this.props.plot.axes[axis] == this.props.fieldId){
        return true
      }
    return false
  }
  render(){
    let icons
    if (this.props.type == 'variable'){
      icons = (
        <span>
          <SVGIcon sprite={invertIcon} active={this.props.invert} onIconClick={()=>{this.props.toggleInvert(this.props.fieldId,this.props.invert)}}/>
        </span>
      )
    }
    else if (this.props.type == 'category'){
      icons = (
        <span>
          <SVGIcon sprite={invertIcon} active={this.props.invert} onIconClick={()=>{this.props.toggleInvert(this.props.fieldId,this.props.invert)}}/>
        </span>
      )
    }
    if (this.props.type == 'selection'){
      icons = (
        <span>
          <SVGIcon sprite={invertIcon} active={this.props.invert} onIconClick={()=>{this.props.toggleInvert(this.props.fieldId,this.props.invert)}}/>
        </span>
      )
    }
    return (
      <div className={styles.header_buttons}>
        {icons}
      </div>
    );
  }
}

export default FieldBoxHeaderFilterButtons
