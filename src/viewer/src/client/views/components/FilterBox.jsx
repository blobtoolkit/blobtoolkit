import React from 'react'
import styles from './Filters.scss';
import plotStyles from './Plot.scss';
import FilterPreview from './FilterPreview'
import FilterBoxHeader from './FilterBoxHeader'
import FilterDisplayRange from './FilterDisplayRange'
import FilterDisplayCategory from './FilterDisplayCategory'
import ExportButton from './ExportButton'

class FilterBox extends React.Component {
  render(){
    let display
    let exportButtons = (
      <span className={plotStyles.download}>
        <ExportButton view={`${this.props.fieldId}_preview`} element={this.props.fieldId+'.preview'} prefix={this.props.datasetId+'.'+this.props.fieldId+'.preview'} format='svg'/>
        <ExportButton view={`${this.props.fieldId}_preview`} element={this.props.fieldId+'.preview'} prefix={this.props.datasetId+'.'+this.props.fieldId+'.preview'} format='png' size={1000}/>
      </span>
    )
    if (this.props.filterType == 'range'){
      display = <FilterDisplayRange {...this.props}/>
    }
    else if (this.props.filterType == 'category'){
      display = <FilterDisplayCategory {...this.props}/>
    }
    return (
      <div className={styles.outer}>

        {display}
        <FilterBoxHeader {...this.props}/>
        {exportButtons}
      </div>
    );
  }
}

export default FilterBox
