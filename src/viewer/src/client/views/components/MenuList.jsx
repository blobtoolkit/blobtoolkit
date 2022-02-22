import React from 'react'
import { connect } from 'react-redux'
import PropTypes from 'prop-types'
import styles from './Layout.scss';
import { Link } from 'react-router-dom'
import Spinner from './Spinner'
import { createSelector } from 'reselect'
import { fetchIdentifiers } from '../reducers/identifiers'
import { getSelectedDatasetMeta } from '../reducers/dataset'
import { getSearchTerm } from '../reducers/location'
import { getSelectSource, getSelectPolygon } from '../reducers/select'
import { getSelectedList, updateSelectedList, getIdentifiersForList, chooseList } from '../reducers/list'
import ListModal from './ListModal';

const ListItem = ({id,list,params,search,onClick,identifiers,active,onChooseList,meta,summaryStats,buscos,selectSource,selectPolygon}) => {
  let css = styles.menu_item
  if (active) css += ' '+styles.active
  let taxon = meta.taxon_name
  let taxid = meta.taxid
  let obj = {
    id,
    datasetId:meta.id,
    taxon,
    taxid,
    search,
    params,
    identifiers
  }
  if (selectSource == 'circle'){
    obj.polygon = selectPolygon
  }
  if (id == 'current'){
    if (!summaryStats.stats){
      return null
    // } ||
    //     !summaryStats.hits || !summaryStats.hits.total ||
    //     Object.keys(summaryStats.busco).length != buscos ||
    //     !summaryStats.taxonomy || !summaryStats.taxonomy.target){
    //   return null
    }
    obj.summaryStats = summaryStats
  }
  return (
    <div id={'list_'+id} className={css}>
      <ListModal name={id} selected={active} dismiss={()=>onClick(null)} list={obj} type={meta.record_type} dataset={meta.id} chooseList={onChooseList}>&nbsp;</ListModal>
      <h3>{id}</h3>
      {list && <span className={styles.menu_subtitle}>{list.length} {meta.record_type}</span>}
    </div>
  )
}

class MenuList extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        identifiers: getIdentifiersForList(state,this.props.id),
        active: getSelectedList(state) == this.props.id,
        meta: getSelectedDatasetMeta(state),
        search: getSearchTerm(state),
        selectSource: getSelectSource(state),
        selectPolygon: getSelectPolygon(state)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        onChooseList: (id,select) => {
          dispatch(chooseList(id,select))
        },
        onClick: (id) => {
          dispatch(updateSelectedList(id))
        }
      }
    }
  }

  render(){
    const ConnectedListItem = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(ListItem)
    return (
      <ConnectedListItem {...this.props}/>
    )
  }
}

export default MenuList
