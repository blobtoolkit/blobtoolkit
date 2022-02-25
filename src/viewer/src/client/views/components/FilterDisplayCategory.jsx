import React from 'react'
import { connect } from 'react-redux'
import styles from './Filters.scss';
import { getCategoryListForFieldId } from '../reducers/preview'

class FilterDisplayCategory extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        return getCategoryListForFieldId(state, props.filterId)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
      }
    }
  }

  render(){
    const ConnectedDisplayCategory = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(DisplayCategory)
    return (
      <ConnectedDisplayCategory filterId={this.props.fieldId} {...this.props}/>
    )
  }
}

class DisplayCategory extends React.Component {
  render() {
    let labels = this.props.bins.map((b,i) => {
      let css = styles.category_label;
      css += b.toggled == true ? ' '+styles.toggled : ''
      return (
        <div key={i}
          className={css}>
          <span className={styles.vertical_text}>
            <span className={styles.highlight}>
              {b.id}
            </span>
          </span>
        </div>
      )
    })
    return (
      <div className={styles.labels_container}>
       {labels}
      </div>
    )
  }
}




export default FilterDisplayCategory
