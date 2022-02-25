import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
import layoutStyles from './Layout.scss'
import { fetchRawData } from '../reducers/field'
import { getBuscoData, getBuscoPaths } from '../reducers/summary'
import { ExportButton } from './ExportButton'
import { format as d3Format } from 'd3-format'
import ReactTable from 'react-table'
import { getDatasetID } from '../reducers/location'
import { addRecords, removeRecords, replaceRecords, setSelectSource } from '../reducers/select'
import { CircleAxis } from './CircleAxis'

const RecordSelector = ({sel,val,records,id,toggleSelect}) => {
  let css = layoutStyles.colored_tab
  if (!sel) return null
  let selected = sel.length / records.length
  let add = false
  if (selected == 1){
    css += ' '+layoutStyles['highlight']
  }
  else if (selected > 0){
    css += ' '+layoutStyles['clear']
    let part_css = layoutStyles['partial_tab']
    add = true
    sel = records
    return (
      <span className={part_css} onClick={()=>toggleSelect(sel,add)}>
      </span>
    )
  }
  else {
    css += ' '+layoutStyles['clear']
    add = true
    sel = records
  }
  return (
    <span className={css} onClick={()=>toggleSelect(sel,add)}/>
  )
}

class BuscoData extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      page:0,
      pageSize:10,
      filter:[],
      sort:[],
      size:[],
      more: false
    }
  }

  componentDidMount() {
    if (!this.props.data){
      this.props.fetchData(this.props.id)
    }
  }

  injectThProps(state, rowInfo, column){
    return {style:{display:'none'}}
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
    let idList = this.state.idList.split(/\W+/)
    let records = []
    idList.forEach(busco=>{
      if (this.props.data.ids.hasOwnProperty(busco)){
        records.push(...this.props.data.ids[busco])
      }
    })
    this.props.replaceRecords(records,true)
  }

  render(){
    if (!this.props.data) return null
    let pct = d3Format(".1%")
    let data = []
    let raw = this.props.data
    // console.log(this.props.data.ids)
    let toggleSelect = this.props.toggleSelect
    data.push({sel:this.props.data.selections.t,records:this.props.data.records.t,cat:'BUSCOs',val:raw.total})
    data.push({sel:this.props.data.selections.c,records:this.props.data.records.c,cat:'Complete',val:raw.scores.c,pct:pct(raw.fractions.c),col:'rgb(51, 160, 44)'})
    data.push({sel:this.props.data.selections.d,records:this.props.data.records.d,cat:'Duplicated',val:raw.scores.d,pct:pct(raw.fractions.d),col:'rgb(32, 100, 27)'})
    data.push({sel:this.props.data.selections.f,records:this.props.data.records.f,cat:'Fragmented',val:raw.scores.f,pct:pct(raw.fractions.f),col:'rgb(163, 226, 127)'})
    data.push({cat:'Missing',val:raw.scores.m,pct:pct(raw.fractions.m),col:'rgb(255, 255, 255)'})
    let columns = [{
      Header:this.props.id,
      columns:
      [
        {
          accessor: 'sel',
          Cell: props => <RecordSelector sel={props.original.sel} val={props.original.val} records={props.original.records} id={[props.original._id]} toggleSelect={toggleSelect}/>,
        width: 40,
          resizable: false
        },
        {
          accessor: 'cat'
        },
        {
          accessor: 'val',
          Cell: props => (<span>{props.original.val}{props.original.col && <span
                              className={layoutStyles.colored_tab}
                              style={{backgroundColor:props.original.col,
                                      float:'left',
                                      border:props.original.cat == 'Missing' ? '0.5px solid #eee' : 'none'}}>
                          </span>}</span>),

        },
        {
          accessor: 'pct'
        }
      ]
    }]
    let csv = data.map(o=>'"'+o.cat+'",'+o.val+','+(o.pct||'')).join('\n')
    let fullList = Object.keys(this.props.data.ids).join(', ')
    let exampleList = Object.keys(this.props.data.ids).slice(0,3).join(', ')
    let exportButtons = (
      <div className={styles.list_download}>
        <ExportButton view={`${this.props.id}_busco`} data={fullList} format='txt' prefix={this.props.datasetId+'.'+this.props.id+'.busco'}/>
      </div>
    )
    return (
      <div style={{display:'flex',flexDirection:'row',minHeight:'14em'}}>
        <div style={{flex:'0 0 50%',display:'flex',alignItems:'center',justifyContent:'flex-end',textAlign:'right',padding:'1em'}}>
          <div style={{flexDirection:'column'}}>
            <ReactTable
                data={data}
                columns={columns}
                showPagination={false}
                defaultPageSize={data.length}
                getTheadThProps={this.injectThProps}
              />
            {fullList && <div style={{maxWidth:'25em'}}>
              <span className={styles.expand}
                    onClick={()=>this.toggleForm()}>{this.state.more ? '...Less' : 'More...'}</span>
              <form className={this.state.more ? '' : layoutStyles.hidden} onSubmit={(e)=>this.handleSubmit(e)}>
                <label>
                  <textarea rows="6"
                            cols="53"
                            value={this.state.idList}
                            placeholder={'Enter a list of BUSCO IDs to select, e.g.:\n'+exampleList}
                            style={{display:'block'}}
                            onChange={(e)=>this.handleChange(e)} />
                </label>
                <input type="submit" value="Select" />
                <div style={{textAlign:'left'}}>
                  {exportButtons}
                  <div>List of BUSCOs in {this.props.datasetId}:</div>
                  <div className={styles.id_list}>{fullList}</div></div>
              </form>
            </div>}
          </div>
        </div>
        <BuscoPlot {...this.props} csv={csv}/>
      </div>
    )
  }
}

class BuscoPlot extends React.Component {
  constructor(props) {
    super(props);
  }

  componentDidMount() {
    if (!this.props.data){
      this.props.fetchData(this.props.id)
    }
  }
  render(){
    let side = 500
    let exportButtons = (
      <div className={styles.download}>
        <ExportButton view={`${this.props.id}_busco`} element={'busco_plot'+this.props.id} prefix={this.props.datasetId+'.'+this.props.id+'.busco'} format='svg'/>
        <ExportButton view={`${this.props.id}_busco`} element={'busco_plot'+this.props.id} prefix={this.props.datasetId+'.'+this.props.id+'.busco'} format='png' size={side}/>
        <ExportButton view={`${this.props.id}_busco`} data={this.props.csv} format='csv' prefix={this.props.datasetId+'.'+this.props.id+'.busco'}/>
      </div>
    )
    let viewbox = '0 0 '+side+' '+side
    if (!this.props.data){
      return null
    }
    return (
      <div style={{flex:'0 0 50%',display:'flex',alignItems:'center',justifyContent:'flex-start',padding:'1em',height:'13em'}}>
        <svg id={'busco_plot'+this.props.id}
          ref={(elem) => { this.svg = elem; }}
          style={{width:'11em',height:'11em',fontSize:'14px'}}
          viewBox={viewbox}
          preserveAspectRatio="xMinYMin">
          <g transform='translate(250,250)'>
            <path d={this.props.paths.c}
                  fill={'rgb(51, 160, 44)'}
                  stroke={'none'}/>
            <path d={this.props.paths.d}
                  fill={'rgb(32, 100, 27)'}
                  stroke={'none'}/>
            <path d={this.props.paths.f}
                  fill={'rgb(163, 226, 127)'}
                  stroke={'none'}/>
            <CircleAxis radius={200}/>
          </g>
        </svg>
        {exportButtons}
      </div>
    )
  }
}

class Busco extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        data: getBuscoData(state,this.props.id),
        paths: getBuscoPaths(state,this.props.id),
        datasetId: getDatasetID(state)
      }
    }
    this.mapDispatchToProps = dispatch => {
      return {
        fetchData: (id) => dispatch(fetchRawData(id)),
        toggleSelect: (id,add) => {
          if (add){
            dispatch(setSelectSource('busco'))
            dispatch(addRecords(id))
          }
          else {
            dispatch(setSelectSource('busco'))
            dispatch(removeRecords(id))
          }
        },
        replaceRecords: arr => {
          dispatch(setSelectSource('busco'))
          dispatch(replaceRecords(arr))
        }
      }
    }

  }

  render(){
    const ConnectedBusco = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(BuscoData)
    return <ConnectedBusco {...this.props}/>
  }
}

export default Busco
