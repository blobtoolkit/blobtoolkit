import React, { Component } from "react";
import { connect } from 'react-redux'
import ReactDOM from "react-dom";
import layoutStyles from './Layout.scss';
import styles from './Customise.scss';
//import 'react-table/react-table.css';
import {
  getDatasetColumns,
  setDatasetColumns,
  columnInfo,
  buscoLineages,
  buscoFields,
  hitTaxa } from '../reducers/datasetTable'


class CustomiseComponent extends Component {
  constructor(props) {
    super(props);
    this.state = {
      buscoLineage: false,
      buscoField: false,
      hitTaxon: false,
      hitField: false
    }
  }

  onToggleField(id,position){
    let list = this.props.columns.slice()
    let index = list.indexOf(id)
    if (index === -1){
      list.splice(position,0,id)
    }
    else {
      list.splice(index, 1)
    }
    this.props.setColumns(list)
  }

  buscoSelector(position){
    let addField
    if (this.state.buscoLineage && this.state.buscoField){
      addField = (<span className={styles.addField}
            onClick={()=>{
              this.onToggleField(`${this.state.buscoLineage}-busco-${this.state.buscoField}`,position);
              this.setState({buscoLineage: false, buscoField: false})
            }}>
            add
      </span>)
    }
    return (
      <span className={styles.field_container}>
        |
        <select onChange={(e)=>{this.setState({buscoLineage:e.target.value})}}
                defaultValue='Lineage'>
          <option disabled={true}>Lineage</option>
          {this.props.buscoLineages.map((lineage,i)=>(
            <option key={i} value={lineage}>
              {lineage}
            </option>
          ))}
        </select>
        <select onChange={(e)=>{this.setState({buscoField:e.target.value})}}
                defaultValue='Value'>
          <option disabled={true} value='Value'>Value</option>
          {Object.keys(buscoFields).map((field,i)=>(
            <option key={i} value={field}>
              {buscoFields[field]}
            </option>
          ))}
        </select>
        {addField}
      </span>

    )
  }



  hitSelector(position){
    const hitFields = {
      count: 'Sequences',
      span: 'Span',
      n50: 'N50'
    }
    let addField
    if (this.state.hitTaxon && this.state.hitField){
      addField = (<span className={styles.addField}
            onClick={()=>{
              this.onToggleField(`${this.state.hitTaxon}-hits-${this.state.hitField}`,position);
              this.setState({hitTaxon: false, hitField: false})
            }}>
            add
      </span>)
    }
    return (
      <span className={styles.field_container}>
        |
        <select onChange={(e)=>{this.setState({hitTaxon:e.target.value})}}
                defaultValue='Taxon'>
          <option disabled={true}>Taxon</option>
          {this.props.hitTaxa.map((taxon,i)=>(
            <option key={i} value={taxon}>
              {taxon}
            </option>
          ))}
        </select>
        <select onChange={(e)=>{this.setState({hitField:e.target.value})}}
                defaultValue='Value'>
          <option disabled={true} value='Value'>Value</option>
          {Object.keys(hitFields).map((field,i)=>(
            <option key={i} value={field}>
              {hitFields[field]}
            </option>
          ))}
        </select>
        {addField}
      </span>

    )
  }


  render() {
    let lastActive = 0
    let groups = Object.keys(this.props.info).map((key,i)=>{
      let group = this.props.info[key]
      let fields = []
      Object.keys(group).forEach((id,j)=>{
        let field = group[id]
        let position = lastActive
        if (typeof field.title !== 'function') {
          fields.push(<span className={styles.field_container} key={j} onClick={()=>this.onToggleField(id,position)}>
            <input type='checkbox'
                   checked={field.active}
                   onChange={()=>this.onToggleField(id,position)}>
            </input>
            <span className={styles.field_title}>{field.title}</span>
          </span>)
          if (field.active){
            lastActive++
          }
        }
      })
      let extra
      if (key == 'BUSCO'){
        extra = this.buscoSelector(lastActive)
      }
      else if (key == 'Assembly statistics'){
        extra = this.hitSelector(lastActive)
      }

      return (
        <div key={i}>
          <span className={styles.columnGroup}>{key}: </span>
          {fields}
          {extra}
        </div>
      )
    })
    return (
      <div className={layoutStyles.listContainer}>
        {groups}
      </div>
    );
  }
}

// getTableProps={(props)=>{
//   return props
//   // props.page
//   // props.pageSize
//   // props.sorted
// }}
class CustomiseDatasetTable extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        columns: getDatasetColumns(state),
        info: columnInfo(state),
        buscoLineages: buscoLineages(state),
        buscoFields,
        hitTaxa: hitTaxa(state)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        setColumns: (arr) => dispatch(setDatasetColumns(arr))
      }
    }
  }

  render(){
    const ConnectedCustomise = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(CustomiseComponent)
    return <ConnectedCustomise {...this.props}/>
  }
}

export default CustomiseDatasetTable
