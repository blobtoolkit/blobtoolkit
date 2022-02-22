import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
import { getStaticBuscoCSV } from '../reducers/summary'
import { ExportButton } from './ExportButton'

class BuscoDownload extends React.Component {
  render(){
    let csv = this.props.buscoData.map(line=>(
      line.map(value=>JSON.stringify(value)).join(',')
    )).join('\n')
    let json = {}
    let keys = this.props.buscoData[0].slice(1)
    this.props.buscoData.slice(1).forEach(row=>{
      let data = {}
      keys.forEach((key,i)=>{
        data[key] = row[(i+1)]
      })
      json[row[0]] = data
    })
    return (
      <div style={{display:'flex',flexDirection:'row',minHeight:'3em'}}>
        <div style={{flex:'0 0 50%',display:'flex',alignItems:'right',justifyContent:'flex-end',textAlign:'left',padding:'1em'}}>
          <div className={styles.download_all} style={{marginTop:0}}>
            export all &nbsp;
            <ExportButton view={'all_busco'}
                          data={csv}
                          format='csv'
                          prefix={this.props.datasetId+'.busco'}/>
            <ExportButton view={'all_busco'}
                          data={json}
                          format='json'
                          prefix={this.props.datasetId+'.busco'}/>
          </div>
          <div></div>
        </div>
      </div>
    )
  }
}

class StaticBuscoData extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        buscoData: getStaticBuscoCSV(state)
      }
    }
  }

  render(){
    const ConnectedBuscoData = connect(
      this.mapStateToProps
    )(BuscoDownload)
    return <ConnectedBuscoData {...this.props}/>
  }
}

export default StaticBuscoData
