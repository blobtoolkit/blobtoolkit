import React from 'react'
import styles from './Filters.scss';
import FilterControlRange from './FilterControlRange'
import FilterControlCategory from './FilterControlCategory'

class FilterBoxHeader extends React.Component {
  render(){
    let control
    if (this.props.filterType == 'range'){
      control = <FilterControlRange {...this.props}/>
    }
    if (this.props.filterType == 'category'){
      control = <FilterControlCategory {...this.props}/>
    }
    return (
      <div className={styles.header}>
        {control}
      </div>
    );
  }
}

export default FilterBoxHeader
