import React from 'react'
import { connect } from 'react-redux'
import styles from './Layout.scss'
import { getDatasetIsFetching,getDatasetIsActive } from '../reducers/repository'
import {
  getTopLevelFields,
  getFieldHierarchy,
  getFieldsByParent
} from '../reducers/field'
import Spinner from './Spinner'
import DatasetApplyFilters from './DatasetApplyFilters'
import DatasetCreateList from './DatasetCreateList'
import Field from './Field'
import FieldSet from './FieldSet'
import MainPlot from './MainPlot'
import ToolTips from './ToolTips'
import MenuDataset from './MenuDataset'
import MenuSwitchView from './MenuSwitchView'
import { getDatasetID } from '../reducers/location'


class FieldMenu extends React.Component {
  mapFields(fields){
    return (
      fields.map((field,i) => {
        let jsx
        if (field.hasRecords){
          jsx = <Field key={i} fieldId={field.id}>{field.id}</Field>
        }
        else if (field.id == 'selection'){
          /*
            TODO:
            Selection Filter goes here
          */
          jsx = <Field key={i} fieldId={field.id}>{field.id}</Field>
        }
        if (field.children){
          return (
            <FieldSet
              key={field.id+'_children'}
              title={field.id}>
              {jsx}
              {this.mapFields(field.children)}
            </FieldSet>
          )
        }
        else {
          return jsx
        }
      })
    )
  }

  render(){
    // console.log(this.props)
    if (!this.props.isActive){
      return (
        <div className={styles.fill_parent}>
          Select a dataset to begin.
        </div>
      )
    }
    if (this.props.isActive == 'loading'){
      return (
        <div className={styles.fill_parent}>
          Loading dataset.
        </div>
      )
    }
    let fields = this.mapFields(this.props.fields)
    return (
      <div className={styles.menu_outer}>
        <MenuSwitchView/>
        <MenuDataset
          key={this.props.datasetId}
          id={this.props.datasetId}
          active={false}
          onDatasetClick={()=>{}}
          onDatasetMount={()=>{}}
        />
        <DatasetApplyFilters />
        {fields}
        <ToolTips set='filterMenu'/>
      </div>
    )
  }
}

const mapStateToProps = state => {
  return {
    isActive: getDatasetIsActive(state),
    topLevelFields: getTopLevelFields(state),
    fields: getFieldHierarchy(state),
    datasetId: getDatasetID(state)
  }
}

const MenuFilterMain = connect(
  mapStateToProps
)(FieldMenu)

export default MenuFilterMain
