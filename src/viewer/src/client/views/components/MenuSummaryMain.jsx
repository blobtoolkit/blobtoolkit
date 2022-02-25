import React from 'react';
import { connect } from 'react-redux'
import styles from './Layout.scss'
import { fetchRawData } from '../reducers/field'
import { getSummary,
  getFullSummary,
  getBuscoSets
} from '../reducers/summary'
import MenuSummary from './MenuSummary'
import MenuDataset from './MenuDataset'
import MenuSwitchView from './MenuSwitchView'
import { getDatasetID, getParsedQueryString, setQueryString } from '../reducers/location'
import { queryToStore } from '../querySync'


class SummaryMenu extends React.Component {
  constructor(props) {
    super(props);
    this.state = {called: -1}
  }

  loadBusco() {
    if (this.props.buscoSets && this.props.fullSummary.buscoScores){
      for (let i = 0; i < this.props.buscoSets.length; i++){
        let set = this.props.buscoSets[i].replace('_busco','')
        if (!this.props.fullSummary.buscoScores[set]){
          if (i > this.state.called){
            this.props.fetchBuscoData(this.props.buscoSets[i])
            this.setState({called:i})
          }
          break
        }
      }
    }
  }

  componentDidMount(){
    this.loadBusco()
  }

  componentDidUpdate(){
    this.loadBusco()
  }

  render(){
    return (
      <div className={styles.menu_outer}>
        <MenuSwitchView/>
        <MenuDataset
          key={this.props.datasetId}
          id={this.props.datasetId}
          active={false}
          onDatasetClick={()=>{}}
          onDatasetMount={()=>{}}
        />

      { this.props.bins ?
          <MenuSummary {...this.props}/> :
          'Select a dataset to view summary statistics'
        }
      </div>
    )
  }

};

const mapStateToProps = state => {
  let summary = getSummary(state) || {}
  let fullSummary = getFullSummary(state) || {}
  let datasetId = getDatasetID(state)
  let buscoSets = getBuscoSets(state)
  let parsed = getParsedQueryString(state)
  return {datasetId,...summary,fullSummary,buscoSets,parsed}
}

const mapDispatchToProps = dispatch => {
  return {
    fetchBuscoData: id => dispatch(fetchRawData(id)),
    onChangeOrder: (values, remove) => {
      dispatch(queryToStore({values,remove}))
    },
    changeQueryParams: (obj) => dispatch(setQueryString(Object.keys(obj).map(k=>`${k}=${obj[k]}`).join('&')))
  }
}

const MenuSummaryMain = connect(
  mapStateToProps,
  mapDispatchToProps
)(SummaryMenu)

export default MenuSummaryMain
