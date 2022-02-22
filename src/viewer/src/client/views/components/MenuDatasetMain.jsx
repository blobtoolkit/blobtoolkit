import {
  fetchMeta,
  getAvailableDatasetIds,
  getRepositoryIsFetching,
  refreshStore,
  setDatasetIsActive,
} from "../reducers/repository";
import {
  getDatasetID,
  getSearchTerm,
  getView,
  toggleHash,
  updatePathname,
} from "../reducers/location";

import MenuDataset from "./MenuDataset";
import MenuSwitchView from "./MenuSwitchView";
import React from "react";
import Search from "./Search";
import Spinner from "./Spinner";
import ToolTips from "./ToolTips";
import { connect } from "react-redux";
import styles from "./Layout.scss";

const DatasetMenu = ({
  searchTerm,
  selectedDataset,
  isFetching,
  datasetIds,
  onDatasetMount,
  onDatasetClick,
}) => {
  return (
    <div className={styles.menu_outer}>
      <MenuSwitchView />
      <Search />
      {isFetching ? <Spinner /> : null}
      {datasetIds.length > 0 ? (
        <span
          className={styles.result_count}
          style={{ marginLeft: "1em" }}
          key="span"
        >
          {datasetIds.length + ' datasets match "' + searchTerm + '"'}
        </span>
      ) : null}
      {datasetIds.map((id) => {
        let active = false;
        if (id == selectedDataset) active = true;
        return (
          <MenuDataset
            key={id}
            id={id}
            active={active}
            onDatasetMount={(id) => onDatasetMount(id)}
            onDatasetClick={(id) => onDatasetClick(id)}
          />
        );
      })}
      <ToolTips set="datasetMenu" />
    </div>
  );
};

const mapStateToProps = (state) => {
  return {
    isFetching: getRepositoryIsFetching(state),
    datasetIds: getAvailableDatasetIds(state),
    selectedDataset: getDatasetID(state),
    searchTerm: getSearchTerm(state),
    view: getView(state),
  };
};

const mapDispatchToProps = (dispatch) => {
  return {
    onDatasetClick: (id, view) => {
      dispatch(refreshStore());
      dispatch(setDatasetIsActive(false));
      if (view) {
        dispatch(updatePathname({ dataset: id, [view]: true }));
      } else {
        dispatch(updatePathname({ dataset: id }));
      }
      dispatch(toggleHash("Filters"));
    },
    onDatasetMount: (id) => dispatch(fetchMeta(id)),
  };
};
const MenuDatasetMenu = connect(
  mapStateToProps,
  mapDispatchToProps
)(DatasetMenu);

export default MenuDatasetMenu;
