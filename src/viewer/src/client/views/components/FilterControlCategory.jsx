import React from 'react'
import { connect } from 'react-redux'
import styles from './Filters.scss';
import { getCategoryListForFieldId } from '../reducers/preview'
import { editFilter } from '../reducers/filter'
import { setDisplayFirst, getMainPlot } from '../reducers/plot'
import { queryToStore } from '../querySync'

function remove(array, element) {
    const index = array.indexOf(element);

    if (index !== -1) {
        array.splice(index, 1);
    }
}

class FilterControlCategory extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        return getCategoryListForFieldId(state, props.filterId)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
      //  onUpdateRange: (id,range) => {return dispatch(editFilter({id:id,range:range}))}
        onToggleBin: (index,keys,filter) => {
          filter.toggled[index] = !filter.toggled[index]
          if (filter.toggled[index]){
            keys = keys.concat(filter.keys)
          }
          else {
            keys.forEach(key=>{
              remove(filter.keys,key)
            })
            keys = filter.keys
          }
          if (keys.length > 0){
            dispatch(queryToStore({
              values:{[this.props.filterId+'--Keys']:keys.join(',')},
              action:'FILTER'
            }))
          }
          else {
            dispatch(queryToStore({remove:[this.props.filterId+'--Keys'],action:'FILTER'}))
          }
          dispatch(editFilter({id:this.props.filterId,toggled:filter.toggled,keys:keys}))
        }
      }
    }
  }

  render(){
    const ConnectedControlCategory = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(ControlCategory)
    return (
      <ConnectedControlCategory filterId={this.props.fieldId} {...this.props}/>
    )
  }
}


class ControlCategory extends React.Component {
  resetInvert(){
    let values = {}
    values[this.props.filterId+'--Inv'] = false
    let remove = [
      this.props.filterId+'--Inv'
    ]
    remove.forEach(key=>{
      delete this.props.parsed[key]
    })
    this.props.onChangeAxisRange(values, remove)
    this.props.changeQueryParams(this.props.parsed)
  }
  resetToggled(){
    let values = {}
    values[this.props.filterId+'--Keys'] = ''
    let remove = [
      this.props.filterId+'--Keys'
    ]
    remove.forEach(key=>{
      delete this.props.parsed[key]
    })
    this.props.onChangeAxisRange(values, remove)
    this.props.changeQueryParams(this.props.parsed)
  }

  render() {
    let bins = this.props.bins.map((b,i) => {
      let css = styles.category_toggle;
      css += b.toggled == true ? ' '+styles.toggled : ''
      return (
        <span key={i}
          className={css}
          rel={this.props.plot.cat}
          style={{backgroundColor:b.color}}
          onClick={()=>this.props.onToggleBin(i,b.keys,this.props.filter)}>
          &nbsp;
        </span>
      )
    })
    let reset
    if (this.props.invert ||
        (this.props.toggled.length > 0 && this.props.toggled.includes(true))){
      reset = (<div className={styles.reset}>
                <span className={styles.reset}
                  onClick={()=>{
                    this.resetInvert()
                    this.resetToggled()
                  }}>
                  reset
                </span>
              </div>)
    }
    return (
      <div className={styles.inside} data-tip data-for='category-toggle'>
        {bins}
        {reset}
      </div>
    )
  }
}


export default FilterControlCategory
