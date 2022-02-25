import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
import layoutStyles from './Layout.scss'
import { getStaticBuscoPaths } from '../reducers/summary'
import { ExportButton } from './ExportButton'
import { format as d3Format } from 'd3-format'
import ReactTable from 'react-table'
import { getDatasetID } from '../reducers/location'
import { CircleAxis } from './CircleAxis'



class BuscoData extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      page:0,
      pageSize:10,
      filter:[],
      sort:[],
      size:[]
    }
  }

  injectThProps(state, rowInfo, column){
    return {style:{display:'none'}}
  }

  render(){
    if (!this.props.data) return null
    let pct = d3Format(".1%")
    let data = []
    let raw = this.props.data
    data.push({cat:'BUSCOs',val:raw.t})
    data.push({cat:'Complete',val:raw.c,pct:pct(raw.c/raw.t),col:'rgb(51, 160, 44)'})
    data.push({cat:'Duplicated',val:raw.d,pct:pct(raw.d/raw.t),col:'rgb(32, 100, 27)'})
    data.push({cat:'Fragmented',val:raw.f,pct:pct(raw.f/raw.t),col:'rgb(163, 226, 127)'})
    data.push({cat:'Missing',val:raw.m,pct:pct(raw.m/raw.t)})
    let columns = [{
      Header:this.props.id,
      columns:
      [
        {
          accessor: 'cat'
        },
        {
          accessor: 'val',
          Cell: props => (<span>{props.original.val}{props.original.col && <span
                              className={layoutStyles.colored_tab}
                              style={{backgroundColor:props.original.col,float:'left'}}>
                          </span>}</span>),

        },
        {
          accessor: 'pct'
        }
      ]
    }]
    let csv = data.map(o=>'"'+o.cat+'",'+o.val+','+(o.pct||'')).join('\n')
    return (
      <div style={{display:'flex',flexDirection:'row',minHeight:'14em'}}>
        <div style={{flex:'0 0 50%',display:'flex',alignItems:'center',justifyContent:'flex-end',textAlign:'right',padding:'1em'}}>
          <ReactTable
              data={data}
              columns={columns}
              showPagination={false}
              defaultPageSize={data.length}
              getTheadThProps={this.injectThProps}
            />
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
    let side = 480
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
      <div style={{flex:'0 0 50%',display:'flex',alignItems:'center',justifyContent:'flex-start',padding:'1em'}}>
        <svg id={'busco_plot'+this.props.id}
          ref={(elem) => { this.svg = elem; }}
          style={{width:'11em',height:'11em',fontSize:'14px'}}
          viewBox={viewbox}
          preserveAspectRatio="xMinYMin">
          <g transform='translate(240,240)'>
            <path d={this.props.paths[this.props.id].c}
                  fill={'rgb(51, 160, 44)'}
                  stroke={'none'}/>
            <path d={this.props.paths[this.props.id].d}
                  fill={'rgb(32, 100, 27)'}
                  stroke={'none'}/>
            <path d={this.props.paths[this.props.id].f}
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

class StaticBusco extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        // paths: getStaticBuscoPaths(state,this.props.id),
        paths: getStaticBuscoPaths(state),
        datasetId: getDatasetID(state)
      }
    }

  }

  render(){
    const ConnectedBusco = connect(
      this.mapStateToProps
    )(BuscoData)
    return <ConnectedBusco {...this.props}/>
  }
}

export default StaticBusco
