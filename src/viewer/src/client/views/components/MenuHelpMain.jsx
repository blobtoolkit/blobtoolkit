import React from 'react';
import { connect } from 'react-redux'
import styles from './Layout.scss'
import { getView } from '../reducers/location'
import FAQs from './FAQs'
import MenuSwitchView from './MenuSwitchView'

const MenuHelp = ({match,datasetId,view}) => {
  switch (view || 'blob') {
    case 'blob':
      view = 'Blob text'
      break
    case 'cumulative':
      view = 'Cumulative text'
      break
    case 'snail':
      view = 'Snail text'
      break
    case 'table':
      view = 'Table text'
      break
    case 'treemap':
      view = 'TreeMap text'
      break
  }
  return (
    <div className={styles.menu_outer}>
      <MenuSwitchView/>
      <FAQs/>
    </div>
  )
  return (
    <span>
      <p>
        This is an alpha release, full documentation will be added once the viewer
        is more feature-complete.
      </p>
      Further information and code are available in the project <a href="https://github.com/blobtoolkit/viewer">Github repository</a>.
      <p>
        {view}
      </p>
    </span>
  )
};

class MenuHelpMain extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        view: getView(state)
      }
    }
  }

  render(){
    const ConnectedHelp = connect(
      this.mapStateToProps
    )(MenuHelp)
    return <ConnectedHelp {...this.props}/>
  }
}

export default MenuHelpMain
