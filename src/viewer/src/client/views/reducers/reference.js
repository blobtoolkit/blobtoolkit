import * as d3 from "d3";

import { addFilter, editFilter, getMetaDataForFilter } from "./filter";
import { createAction, handleAction, handleActions } from "redux-actions";
import {
  getDatasetID,
  getHashString,
  getParsedQueryString,
  getQueryValue,
  getStatic,
  getView,
  setQueryString,
} from "./location";
import {
  getDimensionsbyDimensionId,
  getPreviewDimensions,
  setDimension,
} from "./dimension";

import { apiUrl } from "./api";
import { byIdSelectorCreator } from "./selectorCreators";
import { createSelector } from "reselect";
import deep from "deep-get-set";
import { getCatAxis } from "./plot";
import { getSelectedDatasetMeta } from "./dataset";
import { history } from "./history";
import immutableUpdate from "immutable-update";
import qs from "qs";
import store from "../store";

const requestReferenceValues = createAction("REQUEST_REFERENCE_VALUES");
const receiveReferenceValues = createAction("RECEIVE_REFERENCE_VALUES");
const useStoredReferenceValues = createAction("USE_STORED_REFERENCE_VALUES");
const clearReferenceValues = createAction("CLEAR_REFERENCE_VALUES");

export const referenceValues = handleActions(
  {
    REQUEST_REFERENCE_VALUES: (state, action) => state,
    RECEIVE_REFERENCE_VALUES: (state, action) =>
      immutableUpdate(state, {
        byId: { [action.payload.id]: action.payload.values },
        allIds: [...state.allIds, action.payload.id],
      }),
    USE_STORED_REFERENCE_VALUES: (state, action) => state,
    CLEAR_REFERENCE_VALUES: (state, action) => ({
      byId: {},
      allIds: [],
    }),
  },
  {
    byId: {},
    allIds: [],
  }
);

export const getReferenceValues = (state) => state.referenceValues;

export function resetReferenceValues(fieldId) {
  return (dispatch) => {
    let state = store.getState();
    let params = getParsedQueryString(state);
    let hash = getHashString(state);
    let parsed = {};
    Object.keys(params).forEach((param) => {
      let parts = param.split("--");
      if (parts.length == 3) {
        if (fieldId && parts[1] != fieldId) {
          parsed[param] = params[param];
        }
      }
    });
    let search = qs.stringify(parsed);
    history.push({ search, hash });
    dispatch(setQueryString(search));
    dispatch(clearReferenceValues());
  };
}

const updateQueryString = (params, hash, dispatch) => {
  let search = qs.stringify(params);
  history.push({ search, hash });
  dispatch(setQueryString(search));
};

export const setShowReference = createAction("SET_SHOW_REFERENCE");
export const showReference = handleAction(
  "SET_SHOW_REFERENCE",
  (state, action) => action.payload,
  "false"
);
export const getShowReference = (state) =>
  state.showReference == "true" ? true : false;

export function fetchReferenceValues(fieldId, datasetId, last) {
  return (dispatch) => {
    let state = store.getState();
    let params = getParsedQueryString(state);
    let hash = getHashString(state);
    let referenceValues = getReferenceValues(state);
    let id = `${datasetId}--${fieldId}`;
    let param = `${id}--Active`;
    dispatch(requestReferenceValues(id));
    let values = deep(state, ["referenceValues", "byId", id]);
    if (values) {
      dispatch(useStoredReferenceValues(values));
      if (fieldId != "summary") {
        dispatch(fetchReferenceValues("summary", datasetId));
      }
      if (last) {
        Object.keys(last).forEach((key) => {
          params[key] = true;
        });
        updateQueryString(params, hash, dispatch);
        dispatch(setShowReference("true"));
      }
      return Promise.resolve(useStoredReferenceValues(values));
    }
    let url = `${apiUrl}/field/${datasetId}/${fieldId}`;
    if (fieldId == "summary") {
      url = `${apiUrl}/summary/${datasetId}`;
    }
    return fetch(url)
      .then(
        (response) => response.json()
        // error => console.log('An error occured.', error)
      )
      .then((json) => {
        let values;
        if (fieldId == "summary") {
          dispatch(receiveReferenceValues({ id, values: json.summaryStats }));
        } else {
          if (!json.keys || json.keys.length == 0) {
            values = json.values;
          } else {
            values = json.values.map((v) => json.keys[v]);
          }
          dispatch(receiveReferenceValues({ id, values }));
          dispatch(fetchReferenceValues("summary", datasetId));
          if (last) {
            Object.keys(last).forEach((key) => {
              params[key] = true;
            });
            updateQueryString(params, hash, dispatch);
            dispatch(setShowReference("true"));
          }
        }
      });
  };
}

export const addReferenceFields = () => {
  return (dispatch) => {
    let state = store.getState();
    let params = getParsedQueryString(state);
    let refs = [];
    Object.keys(params).forEach((p) => {
      let parts = p.split("--");
      if (parts.length == 3) {
        if (parts[2] == "Active" && params[p] == "true") {
          refs.push(() =>
            dispatch(fetchReferenceValues(parts[1], parts[0], true))
          );
        }
      }
    });
    if (refs.length > 0) {
      let i = 0;
      let loader = setInterval(() => {
        refs[i]();
        i++;
        if (i == refs.length) {
          clearInterval(loader);
        }
      }, 250);
    }
  };
};

export const referenceReducers = {
  referenceValues,
  showReference,
};
