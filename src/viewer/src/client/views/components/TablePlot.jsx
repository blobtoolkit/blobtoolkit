import React from 'react'
import { connect } from 'react-redux'
import plotStyles from './Plot.scss'
import styles from './Layout.scss'
import './style/node_modules.css'
import {
  getTablePage, setTablePage,
  getTablePageSize, setTablePageSize,
  getTableSortField, setTableSortField,
  getTableSortOrder, setTableSortOrder
} from '../reducers/plotParameters'
import { getTableData,getBinnedColors,getTableDataForPage } from '../reducers/summary'
import { getSelectedDatasetMeta } from '../reducers/dataset'
import { getColorPalette } from '../reducers/color'
import { chooseCurrentRecord, getCurrentRecord } from '../reducers/record'
import { addRecords, removeRecords, recordsByName, setSelectSource } from '../reducers/select'
import CumulativePlotBoundary from './CumulativePlotBoundary'
const saveSvgAsPng = require('save-svg-as-png/lib/saveSvgAsPng.js')
import AxisTitle from './AxisTitle'
import ReactTable from 'react-table'
import ExternalLink from './ExternalLink'
import { fetchIdentifiers } from '../reducers/identifiers'
import DownloadCSV from './TableCSV'
import RecordModal from './RecordModal'
// import { getFilteredList } from '../reducers/filter';
// import { getSelectedRecords } from '../reducers/select'

const CategoryLabel = ({index,keys,id,colors,chooseRecord}) => {
  return (
    <span className={styles.fill_parent}
      onClick={()=>{chooseRecord(id)}}>
      <span
          className={styles.colored_tab}
          style={{backgroundColor:colors[index]}}>
      </span>
      {keys[index]}
    </span>
  )
}

const RecordSelector = ({selected,id,toggleSelect}) => {
  let css = styles.colored_tab
  css += selected ? ' ' + styles['highlight'] : ' ' + styles['clear']
  return (
    <span className={css} onClick={()=>toggleSelect(id,!selected)}/>
  )
}

class CSVWrapper extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      clicked:false
    }
  }

  toggleClicked(clicked){
    if (typeof clicked == 'undefined') clicked = !this.state.clicked
    this.setState({clicked})
  }

  render(){
    let download
    if (this.state.clicked){
      download = <DownloadCSV/>
    }
    return (
      <span className={plotStyles.download}
        onMouseDown={()=>this.toggleClicked(true)}
        onMouseUp={()=>this.toggleClicked(false)}>
        <span className={plotStyles.save_svg}>
          &#8681;csv
          {download}
        </span>
      </span>
    )
  }



}


class Table extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      page:props.page,
      pageSize:props.pageSize,
      filter:[],
      sort:[],
      size:[],
      show:false,
      more:false
    }
    // onPageChange={(pageIndex) => {...}} // Called when the page index is changed by the user
    // onPageSizeChange={(pageSize, pageIndex) => {...}} // Called when the pageSize is changed by the user. The resolve page is also sent to maintain approximate position in the data
    // onSortedChange={(newSorted, column, shiftKey) => {...}} // Called when a sortable column header is clicked with the column itself and if the shiftkey was held. If the column is a pivoted column, `column` will be an array of columns
    // onExpandedChange={(newExpanded, index, event) => {...}} // Called when an expander is clicked. Use this to manage `expanded`
    // onFilteredChange={(column, value) => {...}} // Called when a user enters a value into a filter input field or the value passed to the onFiltersChange handler by the Filter option.
    // onResizedChange={(newResized, event) => {...}} // Called when a user clicks on a resizing component (the right edge of a column header)
  }

  handleChange(e) {
    e.preventDefault()
    this.setState({idList: e.target.value});
  }

  toggleForm() {
    this.setState({more: !this.state.more});
  }

  handleSubmit(e){
    e.preventDefault()
    if (!this.state.idList) return
    let idList = this.state.idList.split(/[\s,;]+/)
    this.props.recordsByName(idList)
  }

  generateColumns(fields,keys,plotCategory,binnedColors){
    let variables = { Header:'Variables', columns: [] }
    let categories = { Header:'Categories', columns: [] }
    Object.keys(fields).forEach(id=>{
      let field = fields[id]
      if (field.type == 'category'){
        let cat = {
          id: field.id,
          Header: () => field.name.replace('_',' '),
          accessor: d => d[id]
        }
        if (id == plotCategory){
          cat.Cell = props => (
            <CategoryLabel index={props.value} id={props.row._id} keys={keys[id]} colors={binnedColors} chooseRecord={this.props.chooseRecord}/>
          )
        }
        else {
          cat.Cell = props => keys[id] ? keys[id][props.value] : 0
        }
        categories.columns.push(cat)
      }
      if (field.type == 'variable'){
        variables.columns.push({
          id: field.id,
          Header: () => field.name.replace('_',' '),
          accessor: d => d[id].toLocaleString(),
          minWidth: 50
        })
      }
    })
    return [variables,categories]
  }
  render(){
    if (!this.props.data) return null
    let keys = this.props.data.keys
    let data = this.props.data.values
    let links = this.props.data.links
    let ids = data.map(o=>o._id)
    let names = data.map(o=>o.id)
    let exampleList = names.slice(0,3).join(', ')
    let toggleSelect = this.props.toggleSelect
    let columns = [{
      Header:'',
      columns:
      [
        {
          accessor: 'selected',
          Cell: props => <RecordSelector selected={props.original.sel} id={[props.original._id]} toggleSelect={toggleSelect}/>,
          width: 30,
          resizable: false,
          sortable: true
        },
        {
          Header: '#',
          accessor: '_id',
          width: 65
        },
        {
          Header: 'ID',
          accessor: 'id',
          show:(data[0] && data[0].id) ? true : false
        },
        {
          Header: 'Links',
          id: 'links',
          accessor: d => (
            links.record.map((link,i)=>(
              <ExternalLink key={i} title={link.title} target='_blank' url={link.func(d)}/>
            ))
          ),
          show:(links.record[0] && links.record[0].title) ? true : false
        }
      ]
    }]
    columns = columns.concat(this.generateColumns(this.props.data.fields,this.props.data.keys,this.props.data.plot.cat.meta.id,this.props.binnedColors))
    let page = this.props.page
    let pages = this.props.data.pages
    let record = this.props.meta.record_type || 'contig'
    let pageSize = this.props.pageSize
    return (
      <div className={plotStyles.fill_parent} style={{display:'block'}}>
        <ReactTable
            data={data}
            page={0}
            pageSize={pageSize}
            columns={columns}
            manual
            getPaginationProps={(props)=>{
              let updated = {...props}
              updated.pages = pages
              updated.page = page
              updated.canNext = updated.page < updated.pages - 1
              updated.canPrevious = updated.page > 0
              return updated
            }}
            onSortedChange={(arr)=>{this.props.updateSort(arr);}}
            onPageChange={(newPage)=>{this.props.changePage(newPage)}}
            onPageSizeChange={(newPageSize, pageIndex)=>{this.props.changePageSize(newPageSize,page,pageSize); return (newPageSize,pageIndex)}}
          />
        <CSVWrapper />
        <RecordModal dismiss={this.props.chooseRecord} currentRecord={this.props.currentRecord.id}/>
        {exampleList && <div style={{maxWidth:'50em'}}>
          <span className={plotStyles.expand}
                onClick={()=>this.toggleForm()}>{this.state.more ? '...Less' : 'More...'}</span>
          <form className={this.state.more ? '' : styles.hidden} onSubmit={(e)=>this.handleSubmit(e)}>
            <label>
              <textarea rows="6"
                        cols="100"
                        value={this.state.idList}
                        placeholder={'Enter a list of '+record+' IDs to select, e.g.:\n'+exampleList}
                        style={{display:'block'}}
                        onChange={(e)=>this.handleChange(e)} />
            </label>
            <input type="submit" value="Select" />
          </form>
        </div>}
      </div>
    )
  }
}

class TablePlot extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      let data = getTableDataForPage(state)
      if (!data) return {}
      let selCount = data.values.filter(o=>o.sel==true).length
      let selectAll = (selCount == data.values.length)
      return {
        data,
        palette: getColorPalette(state),
        page: getTablePage(state),
        pageSize: getTablePageSize(state),
        binnedColors: getBinnedColors(state),
        currentRecord: getCurrentRecord(state),
        meta: getSelectedDatasetMeta(state),
        positions: true,
        selectAll
      }
    }
    this.mapDispatchToProps = dispatch => {
      dispatch(fetchIdentifiers())
      return {
        toggleSelect: (id,add) => {
          if (add){
            dispatch(setSelectSource('table'))
            dispatch(addRecords(id))
          }
          else {
            dispatch(setSelectSource('table'))
            dispatch(removeRecords(id))
          }
        },
        changePage: (page) => {
          dispatch(setTablePage(page))
        },
        changePageSize: (newSize,oldSize,oldPage) => {
          dispatch(setTablePageSize(newSize))
          let newPage = Math.floor(oldSize / newSize * oldPage)
          dispatch(setTablePage(newPage))
        },
        updateSort: (arr) => {
          if (arr.length == 1){
            dispatch(setTableSortField(arr[0].id))
            dispatch(setTableSortOrder(arr[0].desc ? 'desc' : 'asc'))
          }
          dispatch(setTablePage(0))
        },
        chooseRecord: recordId => dispatch(chooseCurrentRecord(recordId)),
        recordsByName: arr => {
          dispatch(setSelectSource('table'))
          dispatch(recordsByName(arr))
        }
      }
    }
  }

  render(){
    const ConnectedTable = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(Table)
    return <ConnectedTable {...this.props}/>
  }
}

export default TablePlot
