import React from 'react'
import { connect } from 'react-redux'
import Header from './Header'
import styles from './Layout.scss'
import MenuDatasetMain from './MenuDatasetMain'
import MenuFilterMain from './MenuFilterMain'
import MenuListsMain from './MenuListsMain'
import MenuDisplayMain from './MenuDisplayMain'
import MenuSummaryMain from './MenuSummaryMain'
import MenuHelpMain from './MenuHelpMain'
import { getTopLevelFields } from '../reducers/field'
import { toggleHash, getHashString, getDatasetID, getStatic } from '../reducers/location'


class HeaderLayoutComponent extends React.Component {
  constructor(props) {
    super(props);
  }

  render(){
    let labels = ['Datasets','Filters','Lists','Settings','Summary','Help']
    let tabs = []
    if (!this.props.datasetId){
      labels = ['Datasets','Help']
    }
    if (this.props.isStatic){
      labels = ['Datasets','Settings','Help']
    }
    let activeTab = this.props.activeTab
    labels.forEach(tab=>{
      tabs.push({id:tab,active:activeTab == tab})
    })
    return (
      <Header tabs={tabs} onTabClick={(tab)=>{this.props.toggleHash(tab)}}/>
    )
  }
}

class LayoutHeader extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        activeTab: getHashString(state),
        datasetId: getDatasetID(state),
        isStatic: getStatic(state)
      }
    },
    this.mapDispatchToProps = dispatch => {
      return {
        toggleHash: value => dispatch(toggleHash(value))
      }
    }
  }

  render(){
    const ConnectedLayoutHeader = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(HeaderLayoutComponent)
    return <ConnectedLayoutHeader {...this.props}/>
  }
}

export default LayoutHeader
