import { addAllFields, editField } from "./field";
import { createAction, handleAction, handleActions } from "redux-actions";
import {
  getDatasetID,
  getHashString,
  getInteractive,
  getQueryString,
  getQueryValue,
  getSearchTerm,
  getStatic,
  getView,
  getViews,
  removeHash,
  updatePathname,
  viewsToPathname,
} from "./location";
import { qsDefault, queryToStore } from "../querySync";

import { addReferenceFields } from "./reference";
import { apiUrl } from "./api";
import { createSelector } from "reselect";
import deep from "deep-get-set";
import { editPlot } from "./plot";
import { fetchIdentifiers } from "./identifiers";
import { filterToList } from "./filter";
import { history } from "./history";
import immutableUpdate from "immutable-update";
import qs from "qs";
import shallow from "shallowequal";
import store from "../store";

const requestRepository = createAction("REQUEST_REPOSITORY");
const receiveRepository = createAction(
  "RECEIVE_REPOSITORY",
  (json) => json,
  () => ({ receivedAt: Date.now() })
);

const invalidateDataset = createAction("INVALIDATE_DATASET");
const requestMeta = createAction("REQUEST_META");
const receiveMeta = createAction(
  "RECEIVE_META",
  (json) => json,
  () => ({ receivedAt: Date.now() })
);
const useStoredMeta = createAction("USE_STORED_META");

const defaultState = () => ({
  isInitialised: false,
  isFetching: false,
  allIds: [],
  byId: {},
});

export const availableDatasets = handleActions(
  {
    REQUEST_REPOSITORY: (state, action) =>
      immutableUpdate(state, {
        isFetching: true,
        isInitialised: true,
        didInvalidate: false,
      }),
    RECEIVE_REPOSITORY: (state, action) => ({
      isFetching: false,
      isInitialised: true,
      byId: action.payload.reduce(
        (obj, item) => ((obj[item.id] = item), obj),
        {}
      ),
      allIds: action.payload.map((item) => item.id),
      lastUpdated: action.meta.receivedAt,
    }),
    INVALIDATE_DATASET: (state, action) => state,
    REQUEST_META: (state, action) => state,
    RECEIVE_META: (state, action) => {
      if (!state.byId[action.payload.id] && !action.payload.json.error) {
        let prefix, latest;
        if (
          action.payload.json.assembly.prefix &&
          action.payload.json.assembly.prefix
        ) {
          prefix = action.payload.json.assembly.prefix;
          latest = Object.keys(state.byId).filter(
            (key) =>
              state.byId[key].prefix === action.payload.json.assembly.prefix
          )[0];
          action.payload.json.latest = state.byId[latest]
            ? state.byId[latest].revision
            : latest;
          action.payload.json.prefix = state.byId[latest]
            ? state.byId[latest].prefix
            : prefix;
        } else {
          action.payload.json.latest = action.payload.json.revision;
        }
      }
      return immutableUpdate(state, {
        byId: { [action.payload.id]: action.payload.json },
      });
    },
  },
  defaultState()
);

export const getRepositoryIsFetching = (state) =>
  state.availableDatasets.isFetching;

export const getRepositoryIsInitialised = (state) =>
  deep(state, "availableDatasets.isInitialised") || false;

export const getAvailableDatasetIds = (state) =>
  deep(state, "availableDatasets.allIds") || [];

export const getAvailableDatasets = (state) =>
  deep(state, "availableDatasets.byId") || [];

export const setApiStatus = createAction("SET_API_STATUS");
export const apiStatus = handleAction(
  "SET_API_STATUS",
  (state, action) => action.payload,
  true
);
export const getApiStatus = (state) => state.apiStatus;

export function fetchRepository(searchTerm) {
  return function (dispatch) {
    dispatch(setDatasetIsActive(false));
    dispatch(requestRepository());
    let state = store.getState();
    let views = getViews(state);
    let datasetId = getDatasetID(state);
    let interactive = getInteractive(state);
    let currentTerm = getSearchTerm(state);
    let setSearch = false;
    let reload = false;
    if (!searchTerm) {
      reload = true;
      searchTerm = currentTerm;
      if (!searchTerm) {
        searchTerm = datasetId;
        if (searchTerm) {
          setSearch = true;
        }
      } else {
        setSearch = true;
      }
    } else if (searchTerm != currentTerm) {
      setSearch = true;
    }
    if (setSearch) {
      dispatch(
        updatePathname({
          search: searchTerm,
          ...([views.primary] && { [views.primary]: true }),
          ...(interactive && { interactive: true }),
        })
      );
    }
    dispatch(refreshStore());

    return fetch(apiUrl + "/search/" + searchTerm)
      .then(
        (response) => response.json(),
        (error) => console.log("An error occured.", error)
      )
      .then((json) => {
        if (datasetId) {
          if (!reload) {
            if (json.length != 1) {
              dispatch(
                updatePathname({}, { dataset: true, [views.primary]: true })
              );
            } else if (json[0].id != datasetId) {
              dispatch(
                updatePathname(
                  { dataset: json[0].id },
                  { [views.primary]: true }
                )
              );
            }
          }
        } else if (json.length == 1) {
          dispatch(
            updatePathname({ dataset: json[0].id }, { [views.primary]: true })
          );
        } else {
          dispatch(
            updatePathname({}, { dataset: true, [views.primary]: true })
          );
        }
        dispatch(receiveRepository(json));
      })
      .catch((err) => dispatch(setApiStatus(false)));
  };
}

function clearUnderscores(obj) {
  let str = JSON.stringify(obj);
  str = str.replace(/"_/g, '"');
  return JSON.parse(str);
}

export function fetchMeta(id) {
  return (dispatch) => {
    dispatch(requestMeta(id));
    let json = deep(store.getState(), ["availableDatasets", "byId", id]);
    if (json && json.records) {
      dispatch(useStoredMeta(json));
      return Promise.resolve(useStoredMeta(json));
    }
    return fetch(`${apiUrl}/dataset/id/${id}`)
      .then((response) => response.json())
      .then((json) => {
        json = clearUnderscores(json);
        json.fields.unshift({
          id: "userDefined",
          children: [{ id: "selection" }],
        });
        dispatch(receiveMeta({ id, json }));
      })
      .catch((err) => {
        dispatch(
          receiveMeta({ id, json: { error: `dataset ${id} not found` } })
        );
      });
  };
}

export const refreshStore = createAction("REFRESH");

window.firstLoad = true;

export const setStaticThreshold = createAction("SET_STATIC_THRESHOLD");
export const chooseStaticThreshold = (staticThreshold) => {
  return function (dispatch) {
    let values = { staticThreshold };
    dispatch(queryToStore({ values }));
  };
};
export const staticThreshold = handleAction(
  "SET_STATIC_THRESHOLD",
  (state, action) => action.payload,
  qs.parse((document.location.search || "").replace(/^\?/, ""))
    .staticThreshold || STATIC_THRESHOLD
);
export const getStaticThreshold = (state) => state.staticThreshold;

export const setNohitThreshold = createAction("SET_NOHIT_THRESHOLD");
export const chooseNohitThreshold = (nohitThreshold) => {
  return function (dispatch) {
    let values = { nohitThreshold };
    dispatch(queryToStore({ values }));
  };
};
export const nohitThreshold = handleAction(
  "SET_NOHIT_THRESHOLD",
  (state, action) => action.payload,
  qs.parse((document.location.search || "").replace(/^\?/, ""))
    .nohitThreshold || NOHIT_THRESHOLD
);
export const getNohitThreshold = (state) => state.nohitThreshold;

export const setMaxSpan = createAction("SET_MAX_SPAN");
export const setMaxCount = createAction("SET_MAX_COUNT");
export const loadDataset = (id, clear) => {
  return function (dispatch) {
    let state = store.getState();
    let isStatic = getStatic(state);
    let interactive = getInteractive(state);
    let threshold = getStaticThreshold(state);
    let nohit = getNohitThreshold(state);
    dispatch(setDatasetIsActive("loading"));
    if (!window.firstLoad && !clear) {
      dispatch(refreshStore());
      let values = {};
      if (threshold != STATIC_THRESHOLD) values.staticThreshold = threshold;
      if (nohit != NOHIT_THRESHOLD) values.nohitThreshold = nohit;
      dispatch(queryToStore({ values, searchReplace: true }));
    }
    dispatch(fetchMeta(id)).then(() => {
      let meta = deep(store.getState(), ["availableDatasets", "byId", id]);
      if (meta.error) {
        console.log(`ERROR: ${meta.error}`);
        dispatch(updatePathname({ notfound: true }, { static: true }));
        dispatch(removeHash());
        dispatch(setDatasetIsActive(true));
      } else {
        dispatch(setMaxCount(meta.assembly["scaffold-count"]));
        dispatch(setMaxSpan(meta.assembly.span));
        let plot = {};
        Object.keys(meta.plot).forEach((key) => {
          plot[key] = meta.plot[key];
        });
        window.plot = plot;
        plot.id = "default";
        Object.keys(plot).forEach((key) => {
          let qv = getQueryValue(key + "Field");
          if (qv) {
            plot[key] = qv;
          }
        });
        let view = getView(state);
        if (
          !interactive &&
          (window.firstLoad || window.records < meta.records) &&
          meta.records > threshold
        ) {
          if (meta.static_plots && !isStatic) {
            dispatch(updatePathname({ [view]: true, static: true }));
          } else if (!meta.static_plots) {
            dispatch(updatePathname({ [view]: true }, { static: true }));
          }
        } else if (!meta.static_plots) {
          dispatch(updatePathname({ [view]: true }, { static: true }));
        } else if (
          window.records > meta.records &&
          meta.records < threshold &&
          isStatic
        ) {
          dispatch(updatePathname({ [view]: true }, { static: true }));
        } else if (meta.records < threshold) {
          dispatch(updatePathname({ [view]: true }, { static: true }));
        }
        dispatch(editPlot(plot));
        Promise.all(addAllFields(dispatch, meta.fields, 1, meta, plot, false))
          .then(() => {
            if (window.firstLoad || clear) {
              let values = qs.parse(getQueryString(state));
              Object.keys(values).forEach((key) => {
                if (key.match("--LimitMax")) {
                  let k = key.replace("--Limit", "--");
                  if (!values.hasOwnProperty(k)) {
                    values[k] = values[key];
                  }
                }
              });
              dispatch(queryToStore({ values }));
              window.firstLoad = false;
              dispatch(filterToList());
            } else {
              dispatch(filterToList());
            }
            window.scrollTop = {};
            window.records = meta.records;
          })
          .then(() => {
            dispatch(setDatasetIsActive(true));
            dispatch(fetchIdentifiers());
            dispatch(addReferenceFields());
          });
      }
    });
  };
};

export const getDatasetMeta = (state, id) =>
  deep(state, ["availableDatasets", "byId", id]) || {};
export const getDatasetIsFetching = (state) => false; //(deep(state,['selectedDataset']) == null) || false

export const setDatasetIsActive = createAction("SET_DATASET_IS_ACTIVE");
export const datasetIsActive = handleAction(
  "SET_DATASET_IS_ACTIVE",
  (state, action) => action.payload,
  false
);

export function fetchSearchResults(str) {
  return function (dispatch) {
    fetch(apiUrl + "/search/" + str).then(
      (response) => response.json(),
      (error) => console.log("An error occured.", error)
    );
  };
}

export const getDatasetIsActive = (state) => {
  return state.datasetIsActive;
};

export const setReloading = createAction("RELOADING");
export const reloading = handleAction(
  "RELOADING",
  (state, action) => action.payload,
  false
);
export const getReloading = (state) => state.reloading;

export const repositoryReducers = {
  availableDatasets,
  fetchRepository,
  apiStatus,
  datasetIsActive,
  reloading,
  staticThreshold,
  nohitThreshold,
};
