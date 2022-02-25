import { createAction, handleAction, handleActions } from "redux-actions";
import { getAnalytics, trackPage } from "./tracking";
import { qsDefault, queryToStore } from "../querySync";
import {
  setDatasetPage,
  setDatasetPageSize,
  setDatasetSorted,
} from "./datasetTable";

import { byIdSelectorCreator } from "./selectorCreators";
import { createSelector } from "reselect";
import { getStaticFields } from "./dataset";
import history from "./history";
import { loadDataset } from "./repository";
import qs from "qs";
import store from "../store";

const basename = BASENAME || "";

export const setPathname = createAction("SET_PATHNAME");
export const pathname = handleAction(
  "SET_PATHNAME",
  (state, action) => action.payload,
  document.location.pathname
    .replace(new RegExp("^" + basename), "")
    .replace("notfound", "") || ""
);
export const getPathname = (state) => {
  return state.pathname;
};

const options = [
  "blob",
  "busco",
  "cumulative",
  "detail",
  "dataset",
  "notfound",
  "report",
  "snail",
  "table",
  "treemap",
  "interactive",
  "static",
];

export const getViews = createSelector(getPathname, (pathname) => {
  let path = pathname.replace(/^\//, "").replace(/\/$/, "").split("/");
  let views = {};
  let nextView = undefined;
  let primary = undefined;
  for (let i = 0; i < path.length; i++) {
    if (options.includes(path[i])) {
      views[path[i]] = true;
      if (
        !primary &&
        path[i] != "search" &&
        path[i] != "dataset" &&
        path[i] != "static"
      ) {
        primary = path[i];
      }
      nextView = path[i];
    } else if (nextView == "dataset") {
      views["dataset"] = path[i];
      nextView = undefined;
    } else if (path[i] !== "undefined") {
      views["search"] = path[i];
    }
  }

  if (views["dataset"] && !views["search"]) {
    views["search"] = views["dataset"];
  }
  views.primary = primary || "blob";
  return views;
});

export const getView = createSelector(getViews, (views) => views.primary);

export const getStatic = createSelector(getViews, (views) => views.static);

export const getInteractive = createSelector(
  getViews,
  (views) => views.interactive
);

export const getDatasetID = createSelector(getViews, (views) => views.dataset);

export const getSearchTerm = createSelector(getViews, (views) =>
  decodeURI(views.search)
);

export const viewsToPathname = (views) => {
  let pathname = "";
  if (views["search"]) {
    pathname += "/" + views["search"];
  }
  if (views["dataset"]) {
    pathname += "/dataset/" + views["dataset"];
  }
  options.forEach((view) => {
    if (view != "dataset" && views[view]) {
      pathname += "/" + view;
    }
  });
  return pathname;
};

export const updatePathname = (update = {}, remove = {}) => {
  return function (dispatch) {
    let state = store.getState();
    let currentPathname = getPathname(state);
    let hash = getHashString(state);
    let search = getQueryString(state);
    let id = getDatasetID(state);
    let views = getViews(state);
    let static_url = state.staticURL;
    let newViews = {};
    if (views.search) newViews.search = views.search;
    if (views.dataset) {
      newViews.dataset = views.dataset;
      if (!newViews.search || newViews.search === "undefined")
        newViews.search = views.dataset;
    }
    let primary = views.primary;
    Object.keys(update).forEach((key) => {
      if (key == "search") {
        dispatch(setDatasetPage(0));
        // dispatch(setDatasetPageSize(10))
        dispatch(setDatasetSorted([]));
      }
      newViews[key] = update[key];
      if (key != "dataset" && key != "static" && options.indexOf(key) != -1) {
        primary = key;
      }
    });
    Object.keys(remove).forEach((key) => {
      delete newViews[key];
      if (key == primary) primary = false;
    });
    if (primary) {
      newViews[primary] = true;
      if (views.static && !remove.static) newViews.static = true;
    }
    let pathname = viewsToPathname(newViews);
    if (pathname != currentPathname) {
      history.push({ pathname, hash, search });
      if (getAnalytics(state)) {
        trackPage(
          pathname + (search ? `?${search}` : "") + (hash ? `#${hash}` : "")
        );
      }
      dispatch(setPathname(pathname));
    }
  };
};

export const chooseView = (view) => {
  return function (dispatch) {
    let state = store.getState();
    let interactive = getInteractive(state);
    dispatch(updatePathname({ [view]: true }));
  };
};

export const toggleStatic = (view, datasetId, isStatic) => {
  return async function (dispatch) {
    if (isStatic) {
      dispatch(updatePathname({ [view]: true }, { static: true }));
      dispatch(loadDataset(datasetId));
    } else {
      dispatch(
        updatePathname({ [view]: true, static: true }, { interactive: true })
      );
    }
  };
};

export const setQueryString = createAction("SET_QUERY_STRING");
export const queryString = handleAction(
  "SET_QUERY_STRING",
  (state, action) => action.payload,
  (document.location.search || "").replace("?", "")
);
export const getQueryString = (state) => state.queryString || "";

export const getParsedQueryString = createSelector(getQueryString, (str) =>
  qs.parse(str)
);

export const setHashString = createAction("SET_HASH_VALUE");
export const hashString = handleAction(
  "SET_HASH_VALUE",
  (state, action) => action.payload,
  (document.location.hash || "").replace("#", "")
);
export const getHashString = (state) => state.hashString || "";

export const toggleHash = (value) => {
  return function (dispatch) {
    let state = store.getState();
    let currentHash = getHashString(state);
    let currentQuery = getQueryString(state);
    if (currentHash && currentHash == value) {
      history.push({ hash: "", search: currentQuery });
      dispatch(setHashString(""));
    } else {
      history.push({ hash: "#" + value, search: currentQuery });
      dispatch(setHashString(value));
    }
  };
};

export const removeHash = (value) => {
  return function (dispatch) {
    let state = store.getState();
    let currentHash = getHashString(state);
    if (value) {
      if (currentHash == value) {
        dispatch(toggleHash(value));
      }
    } else {
      dispatch(toggleHash(currentHash));
    }
  };
};

window.onpopstate = (e) => {
  let state = store.getState();
  let dataset = getDatasetID(state);
  store.dispatch(setPathname(document.location.pathname));
  let newDataset = getDatasetID(state);
  if (newDataset != dataset) {
    store.dispatch(loadDataset(newDataset));
  }
  //store.dispatch(setHashString(document.location.hash.replace(/^#/,'')))
  let currentQuery = getQueryString(state);
  let str = document.location.search.replace(/^\?/, "");
  let values = qs.parse(str);
  store.dispatch(
    queryToStore({ values, searchReplace: true, currentQuery, action: "POP" })
  );
  store.dispatch(setQueryString(str));
};

export const parseQueryString = createSelector(getQueryString, (str = "") => {
  return qs.parse(str.replace("?", ""));
});

const createSelectorForQueryValue = byIdSelectorCreator();

const _getQueryIdAsMemoKey = (queryId) => queryId;

const getQueryId = (queryId) => queryId;

export const getQueryValue = createSelector(
  getQueryId,
  parseQueryString,
  (id, parsed) => {
    return parsed[id] || "";
  }
);

export const locationReducers = {
  pathname,
  hashString,
  queryString,
};
