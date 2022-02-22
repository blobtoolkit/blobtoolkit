import React from 'react'
import { connect } from 'react-redux'
import plotStyles from './Plot.scss'
import styles from './Layout.scss'
import './style/node_modules.css'
import { getSelectedDatasetMeta } from '../reducers/dataset'
import { getSelectedDatasetTable, getLinks } from '../reducers/summary'
import ReactTable from 'react-table'
import ExternalLink from './ExternalLink'
import ExportButton from './ExportButton'
import DownloadCSV from './TableCSV'

class Detail extends React.Component {
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
    let data = this.props.data.data
    let meta = this.props.data.meta
    let columns = [{
      Header:'Assembly metadata',
      columns:
      [
        {
          accessor: 'group',
          width: 100
        },
        {
          accessor: 'key',
          width: 160
        },
        {
          id: 'value',
          accessor: d => {
            let link
            if (d.link){
              let links = Array.isArray(d.link) ? d.link : [d.link]
              link = (
                <span className={styles.detail_link}>
                  {links.map((l,i)=>(
                    <ExternalLink key={i} title={l.title} target='_blank' url={l.func(d.meta)}/>
                  ))}
                </span>
              )
            }
            let value = d.value
            if (d.key == 'span' || d.key == 'scaffold-count'){
              value = value.toLocaleString()
            }
            return (
              <span>
                {link}
                {value}
              </span>
            )
          },
          width: 300
        },
        // {
        //   Header: 'Links',
        //   id: 'links',
        //   accessor: d => (
        //     links.record.map((link,i)=>(
        //       <ExternalLink key={i} title={link.title} target='_blank' url={link.func(d)}/>
        //     ))
        //   ),
        //   show:(links.record[0] && links.record[0].title) ? true : false
        // }
      ]
    }]
    return (
      <div className={plotStyles.outer}
           style={{ top:'2em',
                    maxWidth:'calc(100% - 2em)',
                    maxHeight:'calc(100% - 4em)',
                    overflow:'scroll',
                    display: 'flex',
                    flexDirection: 'column',
                    alignItems: 'center'}}>
        <ReactTable
            data={data}
            page={this.state.page}
            onPageChange={(page)=>{this.setState({page})}}
            columns={columns}
            showPagination={false}
            defaultPageSize={data.length}
            getTheadThProps={this.injectThProps}
            style={{width:'560px'}}
          />
          <span className={plotStyles.download}>
            <ExportButton view='detail' data={meta} format='json' prefix={meta.id+'.meta'}/>
          </span>
      </div>
    )
  }
}

class DetailPlot extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        data: getSelectedDatasetTable(state),
        links: getLinks(state)
      }
    }
  }

  render(){
    const ConnectedDetail = connect(
      this.mapStateToProps
    )(Detail)
    return <ConnectedDetail {...this.props}/>
  }
}

export default DetailPlot
