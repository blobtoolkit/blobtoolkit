import React from 'react'
import { connect } from 'react-redux'
import { getReferenceValues,
  resetReferenceValues,
  fetchReferenceValues,
  getShowReference,
  setShowReference } from '../reducers/reference'
import { getAvailableDatasetIds, getAvailableDatasets } from '../reducers/repository'
import { getSearchTerm } from '../reducers/location'
import { getReferenceLengths } from '../reducers/summary'
import styles from './Layout.scss'
import SVGIcon from './SVGIcon'
import ellipsisIcon from './svg/ellipsis.svg'

export default class ReferenceAssembliesList extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => ({
        current: getReferenceValues(state),
        referenceLengths: getReferenceLengths(state),
        datasetIds: getAvailableDatasetIds(state),
        datasets: getAvailableDatasets(state),
        searchTerm: getSearchTerm(state),
        showReference: getShowReference(state)
      })
    },
    this.mapDispatchToProps = () => {
      return (dispatch) => ({
        reset: (field) => dispatch(resetReferenceValues(field)),
        toggleReference: (value) => dispatch(setShowReference(value)),
        fetchValues: (field, id, last) => dispatch(fetchReferenceValues(field, id, last))
      })
    }
  }

  render(){
    const ConnectedReferenceList = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(ReferenceList)
    return (
      <ConnectedReferenceList {...this.props} />
    )
  }
}

class ReferenceList extends React.Component {
  constructor(props) {
    super(props);
    let ids = Object.keys(this.props.current.byId)
    let expand = ids.length > 0
    ids = (ids.map(id=>id.split('--')[0]) || []).join(', ')
    this.state = {ids, expand}
  }

  handleChange(e) {
    e.preventDefault()
    this.setState({ids: e.target.value});
  }

  toggleForm() {
    this.props.toggleReference(String(!this.state.expand))
    this.setState({expand: !this.state.expand});
  }

  componentDidUpdate(nextProps){
    let ids = Object.keys((this.props.referenceLengths.datasets || {}))
    let newIds = Object.keys((nextProps.referenceLengths.datasets || {}))
    if (ids.join('') != newIds.join('')){
      let expand = newIds.length > 0
      if (expand){

      }
      ids = (newIds.map(id=>id.split('--')[0]) || []).join(', ')
      this.setState({ids, expand})
    }
  }

  handleSubmit(e, list){
    e.preventDefault()
    let ids = []
    if (list){
      ids = list.slice(0)
    }
    else {
        if (!this.state.ids){
        this.props.reset('length')
        this.props.toggleReference('false')
        return
      }
      ids = this.state.ids.trim().split(/[\W+\.-]/)
    }
    this.updateReferenceList(ids, 'length')
    this.props.toggleReference('true')
  }

  updateReferenceList(ids, field){
    this.props.reset(field)
    let list = ids.filter(x=>(x && x.length > 0))
    let params = {}
    list.forEach((id,i)=>{
      let last
      params[`${id}--${field}--Active`] = true
      if (i == list.length - 1){
        last = params
      }
      this.props.fetchValues(field, id, last)
    })
  }

  render(){
    let loadAll
    let list
    let count = this.props.datasetIds.length
    if (count > 1){
      loadAll = `Load all (${this.props.datasetIds.length}) `
      if (this.props.searchTerm != 'all'){
        loadAll += `${this.props.searchTerm} `
      }
      loadAll += 'datasets as reference assemblies'
      list = this.props.datasetIds
    }
    return (
      <div>
        <div className={styles.simple}>
          <h1 className={styles.inline}>reference assemblies</h1>
          <span className={styles.simple_buttons}>
            <SVGIcon sprite={ellipsisIcon} active={this.state.expand} onIconClick={()=>this.toggleForm()}/>
          </span>
        </div>
        <div>
          <form className={this.state.expand ? '' : styles.hidden} onSubmit={(e)=>this.handleSubmit(e)}>
            <input type="submit" value="update" className={styles.menuButton}/>
            <label>
              <textarea rows="3"
                        value={this.state.ids}
                        placeholder={'Enter a list dataset IDs'}
                        className={styles.menuTextarea}
                        onChange={(e)=>this.handleChange(e)} />
            </label>
          </form>
          {loadAll &&<div className={this.state.expand ? '' : styles.hidden}>
            <span onClick={(e)=>this.handleSubmit(e,list)}
                  style={{margin:'0.5em',cursor:'pointer',textDecoration:'underline'}}>
              {loadAll}
            </span>
          </div>}
        </div>
      </div>

    )
  }
}
