import { connect } from 'react-redux'
import { getRepositoryIsFetching, getAvailableDatasetIds } from '../reducers/repository'
import RepositoryDatasetList from '../components/RepositoryDatasetList'
//import { BrowserHistory } from 'react-history'
import { fetchMeta } from '../reducers/repository'
//import Spinner from '../components/Spinner'

const mapStateToProps = state => {
  return {
    isFetching: getRepositoryIsFetching(state),
    datasetIds: getAvailableDatasetIds(state)
  }
}

const mapDispatchToProps = dispatch => {
  return {
    // onDatasetClick: id => dispatch(selectedDataset(id)),
    onDatasetMount: id => dispatch(fetchMeta(id))
  }
}
const Repository = connect(
  mapStateToProps,
  mapDispatchToProps
)(RepositoryDatasetList)

export default Repository
