import { createAction, handleAction, handleActions } from "redux-actions";
import { defaultTransform, qsDefault, queryToStore } from "../querySync";
import { fetchRawData, getRawDataForFieldId } from "./field";
import { pathname, setPathname, urlViews, viewsToPathname } from "./repository";

import { byIdSelectorCreator } from "./selectorCreators";
import { createSelector } from "reselect";
import { getDatasetID } from "./location";
import { getFilteredList } from "./filter";
import { getMainPlot } from "./plot";
import history from "./history";
import immutableUpdate from "immutable-update";
import qs from "qs";
import store from "../store";

export const setPlotShape = createAction("SET_PLOT_SHAPE");
export const choosePlotShape = (plotShape) => {
  return function (dispatch) {
    let values = { plotShape };
    dispatch(queryToStore({ values }));
  };
};
export const plotShape = handleAction(
  "SET_PLOT_SHAPE",
  (state, action) => action.payload,
  qsDefault("plotShape")
);
export const getPlotShape = (state) => state.plotShape;

export const setPlotResolution = createAction("SET_PLOT_RESOLUTION");
export const choosePlotResolution = (plotResolution) => {
  return function (dispatch) {
    let values = { plotResolution };
    dispatch(queryToStore({ values }));
  };
};
export const plotResolution = handleAction(
  "SET_PLOT_RESOLUTION",
  (state, action) => action.payload,
  qsDefault("plotResolution")
);
export const getPlotResolution = (state) => state.plotResolution;

export const setPlotScale = createAction("SET_PLOT_SCALE");
export const choosePlotScale = (plotScale) => {
  return function (dispatch) {
    let values = { plotScale };
    dispatch(queryToStore({ values }));
  };
};
export const plotScale = handleAction(
  "SET_PLOT_SCALE",
  (state, action) => action.payload,
  qsDefault("plotScale")
);
export const getPlotScale = (state) => state.plotScale;

export const setZScale = createAction("SET_Z_SCALE");
export const chooseZScale = (zScale) => {
  return function (dispatch) {
    let values = { zScale };
    dispatch(queryToStore({ values }));
  };
};
export const zScale = handleAction(
  "SET_Z_SCALE",
  (state, action) => action.payload,
  qsDefault("zScale")
);
export const getZScale = (state) => state.zScale;

export const setZReducer = createAction("SET_Z_REDUCER");
export const chooseZReducer = (zReducer) => {
  return function (dispatch) {
    let values = { zReducer };
    dispatch(queryToStore({ values }));
  };
};
export const zReducers = {
  sum: (arr) => {
    let sum = 0;
    let len = arr.length;
    for (let i = 0; i < len; i++) {
      sum += arr[i];
    }
    return sum;
  },
  min: (arr) => {
    let min = Number.POSITIVE_INFINITY;
    let len = arr.length;
    for (let i = 0; i < len; i++) {
      if (arr[i] < min) min = arr[i];
    }
    return min;
  },
  max: (arr) => {
    let max = Number.NEGATIVE_INFINITY;
    let len = arr.length;
    for (let i = 0; i < len; i++) {
      if (arr[i] > max) max = arr[i];
    }
    return max;
  },
  count: (arr) => (arr.length > 0 ? arr.length : 0),
  mean: (arr) => {
    if (!arr.length > 0) return 0;
    let sum = 0;
    let len = arr.length;
    for (let i = 0; i < len; i++) {
      sum += arr[i];
    }
    return sum / arr.length;
  },
  n50: (arr) => {
    if (!arr.length > 0) return 0;
    arr.sort((a, b) => b - a);
    let len = arr.length;
    let sum = 0;
    for (let i = 0; i < len; i++) {
      sum += arr[i];
    }
    let csum = 0;
    for (let i = 0; i < len; i++) {
      csum += arr[i];
      if (csum >= Math.floor(sum / 2)) {
        return arr[i];
      }
    }
    return 0;
  },
};
export const zReducer = handleAction(
  "SET_Z_REDUCER",
  (state, action) => ({ id: action.payload, func: zReducers[action.payload] }),
  { id: "sum", func: zReducers.sum }
);
export const getZReducer = (state) => state.zReducer;

export const setTransformFunction = createAction("EDIT_Y_TRANSFORM");
export const transformFunction = handleAction(
  "EDIT_Y_TRANSFORM",
  (state, action) => {
    let obj = {};
    obj.intercept = action.payload.hasOwnProperty("intercept")
      ? action.payload.intercept * 1
      : defaultTransform.intercept;
    obj.x = action.payload.hasOwnProperty("x")
      ? action.payload.x * 1
      : defaultTransform.x;
    obj.order = action.payload.hasOwnProperty("order")
      ? action.payload.order * 1
      : defaultTransform.order;
    obj.factor = action.payload.hasOwnProperty("factor")
      ? action.payload.factor * 1
      : defaultTransform.factor;
    obj.index = action.payload.hasOwnProperty("index")
      ? action.payload.index * 1
      : defaultTransform.index;
    obj.origin = action.payload.hasOwnProperty("origin")
      ? action.payload.origin
      : defaultTransform.origin;
    return obj;
  },
  defaultTransform
);
export const getTransformFunctionParams = (state) => state.transformFunction;

export const getTransformFunction = createSelector(
  getTransformFunctionParams,
  (props) => {
    let factor = props.factor / 900;
    let intercept = props.intercept;
    return ([x, y]) => {
      let newY =
        y +
        Math.abs(props.x - x) ** props.order *
          (props.factor / 900 ** (props.order - 1)) +
        props.intercept;
      return [x, newY];
    };
  }
);
export const chooseTransformFunction = (obj) => {
  return function (dispatch) {
    let values = {};
    let remove = [];
    let origin = { x: 0, y: 0 };
    Object.keys(obj).forEach((key) => {
      if (key == "intercept" || key == "factor") {
        let k = `transform${key[0].toUpperCase() + key.slice(1)}`;
        if (obj[key] != defaultTransform[key] && !isNaN(obj[key])) {
          values[k] = obj[key];
        } else {
          remove.push(k);
        }
      }
    });
    if (Object.keys(obj).length == 0) {
      Object.keys(defaultTransform).forEach((key) => {
        if (key == "intercept" || key == "factor") {
          let k = `transform${key[0].toUpperCase() + key.slice(1)}`;
          values[k] = defaultTransform[key];
          remove.push(k);
        }
      });
    }
    dispatch(queryToStore({ values, remove }));
  };
};

export const setCurveOrigin = createAction("SET_CURVE_ORIGIN");
export const chooseCurveOrigin = (curveOrigin) => {
  return function (dispatch) {
    let values = { curveOrigin };
    dispatch(queryToStore({ values }));
  };
};
export const curveOrigin = handleAction(
  "SET_CURVE_ORIGIN",
  (state, action) => action.payload,
  qsDefault("curveOrigin")
);
export const getCurveOrigin = (state) => state.curveOrigin;

export const setScaleTo = createAction("SET_SCALE_TO");
export const chooseScaleTo = (scaleTo) => {
  return function (dispatch) {
    let values = { scaleTo };
    dispatch(queryToStore({ values }));
  };
};
export const scaleTo = handleAction(
  "SET_SCALE_TO",
  (state, action) => action.payload,
  qsDefault("scaleTo")
);
export const getScaleTo = (state) => state.scaleTo;

export const setCircumferenceScale = createAction("SET_CIRCUMFERENCE_SCALE");
export const chooseCircumferenceScale = (circumferenceScale) => {
  return function (dispatch) {
    let values = { circumferenceScale };
    dispatch(queryToStore({ values }));
  };
};
export const circumferenceScale = handleAction(
  "SET_CIRCUMFERENCE_SCALE",
  (state, action) => action.payload,
  qsDefault("circumferenceScale")
);
export const getCircumferenceScale = (state) => state.circumferenceScale;

export const setRadiusScale = createAction("SET_RADIUS_SCALE");
export const chooseRadiusScale = (radiusScale) => {
  return function (dispatch) {
    let values = { radiusScale };
    dispatch(queryToStore({ values }));
  };
};
export const radiusScale = handleAction(
  "SET_RADIUS_SCALE",
  (state, action) => action.payload,
  qsDefault("radiusScale")
);
export const getRadiusScale = (state) => state.radiusScale;

export const setSnailOrigin = createAction("SET_SNAIL_ORIGIN");
export const chooseSnailOrigin = (snailOrigin) => {
  return function (dispatch) {
    let values = { snailOrigin };
    dispatch(queryToStore({ values }));
  };
};
export const snailOrigin = handleAction(
  "SET_SNAIL_ORIGIN",
  (state, action) => action.payload,
  qsDefault("snailOrigin")
);
export const getSnailOrigin = (state) => state.snailOrigin;

export const setSideMax = createAction("SET_SIDE_SCALE");
export const chooseSideMax = (sideMax) => {
  return function (dispatch) {
    let values = { sideMax };
    dispatch(queryToStore({ values }));
  };
};
export const sideMax = handleAction(
  "SET_SIDE_SCALE",
  (state, action) => action.payload,
  qsDefault("sideMax")
);
export const getSideMax = (state) => state.sideMax;

export const chooseMaxSpan = (maxSpan) => {
  return function (dispatch) {
    let values = { maxSpan };
    dispatch(queryToStore({ values }));
  };
};
export const maxSpan = handleAction(
  "SET_MAX_SPAN",
  (state, action) => action.payload,
  qsDefault("maxSpan")
);
export const getMaxSpan = (state) => state.maxSpan;

export const chooseMaxCount = (maxCount) => {
  return function (dispatch) {
    let values = { maxCount };
    dispatch(queryToStore({ values }));
  };
};
export const maxCount = handleAction(
  "SET_MAX_COUNT",
  (state, action) => action.payload,
  qsDefault("maxCount")
);
export const getMaxCount = (state) => state.maxCount;

export const setSVGThreshold = createAction("SET_SVG_THRESHOLD");
export const chooseSVGThreshold = (svgThreshold) => {
  return function (dispatch) {
    let values = { svgThreshold };
    dispatch(queryToStore({ values }));
  };
};
export const svgThreshold = handleAction(
  "SET_SVG_THRESHOLD",
  (state, action) => action.payload,
  qsDefault("svgThreshold")
);
export const getSVGThreshold = (state) => state.svgThreshold;

export const setPlotGraphics = createAction("SET_PLOT_GRAPHICS");
export const choosePlotGraphics = (plotGraphics) => {
  return function (dispatch) {
    let values = { plotGraphics };
    dispatch(queryToStore({ values }));
  };
};
export const plotGraphics = handleAction(
  "SET_PLOT_GRAPHICS",
  (state, action) => action.payload,
  qsDefault("plotGraphics")
);
const _getPlotGraphics = (state) => state.plotGraphics;

export const getPlotGraphics = createSelector(
  _getPlotGraphics,
  getSVGThreshold,
  getFilteredList,
  (graphics, threshold, list) => {
    if (graphics != "svg" && graphics != "canvas") {
      graphics = "canvas";
      if (threshold >= list.length) {
        graphics = "svg";
      }
    }
    return graphics;
  }
);

export const setCircleLimit = createAction("SET_CIRCLE_LIMIT");
export const chooseCircleLimit = (circleLimit) => {
  return function (dispatch) {
    let values = { circleLimit };
    dispatch(queryToStore({ values }));
  };
};
export const circleLimit = handleAction(
  "SET_CIRCLE_LIMIT",
  (state, action) => action.payload,
  qs.parse((document.location.search || "").replace(/^\?/, "")).circleLimit ||
    CIRCLE_LIMIT
);
export const getCircleLimit = (state) => state.circleLimit;

export const setTablePage = createAction("SET_TABLE_PAGE");
export const tablePage = handleAction(
  "SET_TABLE_PAGE",
  (state, action) => action.payload,
  0
);
export const getTablePage = (state) => state.tablePage;

export const setTablePageSize = createAction("SET_TABLE_PAGE_SIZE");
export const tablePageSize = handleAction(
  "SET_TABLE_PAGE_SIZE",
  (state, action) => action.payload,
  10
);
export const getTablePageSize = (state) => state.tablePageSize;

export const setTableSortField = createAction("SET_TABLE_SORT_FIELD");
export const tableSortField = handleAction(
  "SET_TABLE_SORT_FIELD",
  (state, action) => action.payload,
  null
);
export const getTableSortField = (state) => state.tableSortField;

export const setTableSortOrder = createAction("SET_TABLE_SORT_ORDER");
export const tableSortOrder = handleAction(
  "SET_TABLE_SORT_ORDER",
  (state, action) => action.payload,
  "desc"
);
export const getTableSortOrder = (state) => state.tableSortOrder;

export const setPngResolution = createAction("SET_PNG_RESOLUTION");
export const choosePngResolution = (pngResolution) => {
  return function (dispatch) {
    let values = { pngResolution };
    dispatch(queryToStore({ values }));
  };
};
export const pngResolution = handleAction(
  "SET_PNG_RESOLUTION",
  (state, action) => action.payload,
  qsDefault("pngResolution")
);
export const getPngResolution = (state) => state.pngResolution;

export const setOtherLimit = createAction("SET_OTHER_LIMIT");
export const chooseOtherLimit = (otherLimit) => {
  return function (dispatch) {
    let values = { otherLimit };
    dispatch(queryToStore({ values }));
  };
};
export const otherLimit = handleAction(
  "SET_OTHER_LIMIT",
  (state, action) => action.payload,
  qsDefault("otherLimit")
);
export const getOtherLimit = (state) => state.otherLimit;

export const setShowTotal = createAction("SET_SHOW_TOTAL");
export const chooseShowTotal = (showTotal) => {
  return function (dispatch) {
    showTotal = String(showTotal) == "true" ? "true" : "false";
    let values = { showTotal };
    dispatch(queryToStore({ values }));
  };
};
export const showTotal = handleAction(
  "SET_SHOW_TOTAL",
  (state, action) => action.payload,
  qsDefault("showTotal")
);
export const getShowTotal = (state) =>
  state.showTotal == "true" ? true : false;

export const setLargeFonts = createAction("SET_LARGE_FONTS");
export const chooseLargeFonts = (largeFonts) => {
  return function (dispatch) {
    largeFonts = String(largeFonts) == "true" ? "true" : "false";
    let values = { largeFonts };
    dispatch(queryToStore({ values }));
  };
};
export const largeFonts = handleAction(
  "SET_LARGE_FONTS",
  (state, action) => action.payload,
  qsDefault("largeFonts")
);
export const getLargeFonts = (state) =>
  state.largeFonts == "true" ? true : false;

export const setErrorBars = createAction("SET_ERROR_BARS");
export const chooseErrorBars = (errorBars) => {
  return function (dispatch) {
    let values = { errorBars };
    dispatch(queryToStore({ values }));
  };
};
export const errorBars = handleAction(
  "SET_ERROR_BARS",
  (state, action) => action.payload,
  qsDefault("errorBars")
);
export const getErrorBars = (state) => state.errorBars;

export const setWindowSize = createAction("SET_WINDOW_SIZE");
export const chooseWindowSize = (windowSize) => {
  return function (dispatch) {
    let values = { windowSize };
    dispatch(queryToStore({ values }));
  };
};
export const windowSize = handleAction(
  "SET_WINDOW_SIZE",
  (state, action) => action.payload,
  qsDefault("windowSize")
);
export const getWindowSize = (state) => state.windowSize;

// export const setAdjustCoverage = createAction('SET_ADJUST_COVERAGE')
// export const chooseAdjustCoverage = (adjust) => {
//   return function (dispatch) {
//     adjust = String(adjust) == 'true' ? 'true' : 'false'
//     let values = {adjustCoverage: adjust}
//     if (adjust == 'true'){
//       let state = store.getState()
//       let ncount = getRawDataForFieldId(state, 'ncount')
//       if (!ncount){
//         dispatch(fetchRawData('ncount'))
//       }
//       values['ncount--Active'] = true
//     }
//     dispatch(queryToStore({values}))
//   }
// }
// export const adjustCoverage = handleAction(
//   'SET_ADJUST_COVERAGE',
//   (state, action) => (
//     action.payload
//   ),
//   qsDefault('adjustCoverage')
// )
// export const getAdjustCoverage = state => state.adjustCoverage == 'true' ? true : false

export const plotParameterReducers = {
  plotShape,
  plotResolution,
  plotGraphics,
  svgThreshold,
  pngResolution,
  plotScale,
  zScale,
  zReducer,
  transformFunction,
  curveOrigin,
  scaleTo,
  circumferenceScale,
  radiusScale,
  snailOrigin,
  sideMax,
  maxCount,
  maxSpan,
  circleLimit,
  tablePage,
  tablePageSize,
  tableSortField,
  tableSortOrder,
  otherLimit,
  showTotal,
  // adjustCoverage,
  largeFonts,
  errorBars,
  windowSize,
};
