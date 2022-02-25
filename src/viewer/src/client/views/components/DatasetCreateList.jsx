import React from 'react'
import { connect } from 'react-redux'
import {
  getFilteredList
} from '../reducers/filter'
import {
  createList
} from '../reducers/list'
import styles from './Datasets.scss'


const mapStateToProps = state => ({list:getFilteredList(state)})

const mapDispatchToProps = dispatch => {
  return {
    onClick: (obj) => dispatch(createList(obj))
  }
}

class CreateList extends React.Component {
  constructor(props) {
    super(props);
    this.state = {id: 'default'};
    this.handleChange = this.handleChange.bind(this)
  }
  handleChange(event){
    this.setState({id: event.target.value})
  }
  render(){
    return (
      <div className={styles.dataset_controls}>
        <div className={styles.detail_container + ' ' + styles.list_detail}>
          <input type='text'
            onChange={this.handleChange}
            value={this.state.id}
            data-tip data-for='create-list-input'/>
          <span onClick={()=>this.props.onClick({id:this.state.id,list:this.props.list})}
            className={styles.button + ' ' + styles.list_button}
            data-tip data-for='create-list-button' >
            Create List
          </span>
        </div>
      </div>
    )
  }
}

const DatasetCreateList = connect(
  mapStateToProps,
  mapDispatchToProps
)(CreateList)

export default DatasetCreateList
