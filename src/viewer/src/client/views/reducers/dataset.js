import { createAction, handleAction, handleActions } from "redux-actions";
import { createSelector, createSelectorCreator } from "reselect";

import { byIdSelectorCreator } from "./selectorCreators";
import deep from "deep-get-set";
import { getDatasetID } from "./location";
import immutableUpdate from "immutable-update";
import shallow from "shallowequal";
import store from "../store";

export const getMetaDataForDataset = (state, datasetId) =>
  state.availableDatasets.byId[datasetId] || {};
const _getDatasetIdAsMemoKey = (state, datasetId) => datasetId;

const createSelectorForDatasetMeta = byIdSelectorCreator();

export const getDatasetMetaById = createSelectorForDatasetMeta(
  _getDatasetIdAsMemoKey,
  getMetaDataForDataset,
  (meta) => {
    return meta || {};
  }
);

export const getCurrentDatasetMeta = createSelector(
  (state) => getDatasetMetaById(state, getDatasetID(state)),
  (meta) => meta
);

export const getSelectedDatasetMeta = createSelector(
  getCurrentDatasetMeta,
  (meta) => meta
);

export const getDatasetName = createSelector(getCurrentDatasetMeta, (meta) =>
  meta ? meta.name : ""
);

export const getStaticFields = createSelector(
  getCurrentDatasetMeta,
  (meta) => meta.static_plots || false
);

export const getDatasetsById = (state, datasetId) =>
  state.datasets
    ? state.datasets.byId
      ? state.datasets.byId[datasetId]
      : {}
    : {};

const createSelectorForDatasetId = byIdSelectorCreator();

export const getDatasetById = createSelectorForDatasetId(
  _getDatasetIdAsMemoKey,
  getDatasetsById,
  (dataset) => dataset || {}
);

export const getCurrentDataset = createSelector(
  (state) => getDatasetById(state, getDatasetID(state)),
  (dataset) => dataset
);

export const getDataset = createSelector(
  getCurrentDataset,
  (dataset) => dataset
);
