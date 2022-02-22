import React from 'react'
import PropTypes from 'prop-types'
import styles from './Repository.scss';
import RepositoryDataset from './RepositoryDataset'
import Spinner from './Spinner'

const RepositoryDatasetList = ({ isFetching, datasetIds, onDatasetMount }) => (
  <div>
    {isFetching ? <Spinner /> : null}
    {datasetIds.map(id => (
      <RepositoryDataset key={id} id={id} onMount={(id) => onDatasetMount(id)} />
    ))}
  </div>
)

RepositoryDatasetList.propTypes = {
  isFetching: PropTypes.bool.isRequired,
  datasetIds: PropTypes.array.isRequired
}

export default RepositoryDatasetList
