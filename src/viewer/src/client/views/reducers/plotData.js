import * as d3 from "d3";

import { createAction, handleAction, handleActions } from "redux-actions";
import {
  getBinsForCat,
  getBinsForFieldId,
  getDetailsForFieldId,
  getFields,
  getMetaDataForField,
  getOtherLimit,
  getRawDataForFieldId,
} from "./field";
import { getCatAxis, getMainPlot, getXAxis, getYAxis, getZAxis } from "./plot";
import {
  getCategoryListForFieldId,
  getFilteredDataForFieldId,
  getPlainCategoryListForFieldId,
} from "./preview";
import {
  getErrorBars,
  getPlotResolution,
  getPlotScale,
  getTransformFunction,
  getTransformFunctionParams,
  getWindowSize,
  getZScale,
  radiusScale,
} from "./plotParameters";

import { byIdSelectorCreator } from "./selectorCreators";
import { clampDomain } from "./field";
import cloneFunction from "clone-function";
import { createSelector } from "reselect";
import { datasetColumns } from "./datasetTable";
import { getColorPalette } from "./color";
import { getDatasetID } from "./location";
import { getFilteredList } from "./filter";
import { getQueryValue } from "./location";
import store from "../store";

const getAxis = (state, axis) => axis;

export const getAxisTitle = createSelector(
  getAxis,
  getMainPlot,
  (axis, plot) => {
    return plot.axes[axis];
  }
);

const getRawDataForX = createSelector(
  (state) => getRawDataForFieldId(state, getXAxis(state)),
  (data) => data
);

const getRawDataForY = createSelector(
  (state) => getRawDataForFieldId(state, getYAxis(state)),
  (data) => data
);

const getRawDataForZ = createSelector(
  (state) => getRawDataForFieldId(state, getZAxis(state)),
  (data) => data
);

export const getRawDataForCat = createSelector(
  (state) => getRawDataForFieldId(state, getCatAxis(state)),
  (data) => data
);

const getFilteredDataForX = createSelector(
  (state) => getFilteredDataForFieldId(state, getXAxis(state)),
  (data) => {
    return data;
  }
);

const getFilteredDataForY = createSelector(
  (state) => getFilteredDataForFieldId(state, getYAxis(state)),
  (data) => data
);

const getFilteredDataForZ = createSelector(
  (state) => getFilteredDataForFieldId(state, getZAxis(state)),
  (data) => data
);

export const getFilteredDataForCat = createSelector(
  (state) => getFilteredDataForFieldId(state, getCatAxis(state)),
  (data) => data
);

export const getFilteredDataForGC = createSelector(
  (state) => getFilteredDataForFieldId(state, "gc"),
  (data) => data
);

export const getFilteredDataForLength = createSelector(
  (state) => getFilteredDataForFieldId(state, "length"),
  (data) => data
);

export const getFilteredDataForNCount = createSelector(
  (state) => getFilteredDataForFieldId(state, "ncount"),
  (data) => data
);

export const getDetailsForX = createSelector(
  (state) => getDetailsForFieldId(state, getXAxis(state)),
  (data) => data
);

export const getDetailsForY = createSelector(
  (state) => getDetailsForFieldId(state, getYAxis(state)),
  (data) => data
);

const getDetailsForZ = createSelector(
  (state) => getDetailsForFieldId(state, getZAxis(state)),
  (data) => data
);

export const getDetailsForCat = createSelector(
  (state) => getDetailsForFieldId(state, getCatAxis(state)),
  (data) => data
);

export const getMainPlotData = createSelector(
  getDatasetID,
  getFilteredDataForX,
  getFilteredDataForY,
  getFilteredDataForZ,
  getFilteredDataForCat,
  getDetailsForX,
  getDetailsForY,
  getDetailsForZ,
  getDetailsForCat,
  (dsId, xData, yData, zData, catData, xMeta, yMeta, zMeta, catMeta) => {
    if (!dsId || !catData) return undefined;
    let plotData = { id: "default", axes: {}, meta: {} };
    plotData.axes.x = xData || { values: [] };
    let xDomain = xMeta.xScale ? xMeta.xScale.domain().slice(0) : [0, 10];
    let xmin = getQueryValue("xmin");
    if (xmin) {
      xDomain[0] = 1 * xmin;
    }
    let xmax = getQueryValue("xmax");
    if (xmax) {
      xDomain[1] = 1 * xmax;
    }
    xMeta.xScale.domain(xDomain);
    xMeta.xScale.type = xMeta.meta.scale;
    plotData.meta.x = xMeta;
    plotData.axes.y = yData || { values: [] };
    let yDomain = yMeta.xScale ? yMeta.xScale.domain().slice(0) : [0, 10];
    let ymin = getQueryValue("ymin");
    if (ymin) {
      yDomain[0] = 1 * ymin;
    }
    let ymax = getQueryValue("ymax");
    if (ymax) {
      yDomain[1] = 1 * ymax;
    }
    if (yMeta.xScale) {
      yMeta.xScale.domain(yDomain);
      yMeta.xScale.type = yMeta.meta.scale;
    }
    plotData.meta.y = yMeta;
    plotData.axes.z = zData || { values: [] };
    plotData.meta.z = zMeta;
    plotData.axes.cat = catData;
    plotData.meta.cat = catMeta;
    return plotData;
  }
);

export const getAllMainPlotData = createSelector(
  getDatasetID,
  getRawDataForX,
  getRawDataForY,
  getRawDataForZ,
  getRawDataForCat,
  getDetailsForX,
  getDetailsForY,
  getDetailsForZ,
  getDetailsForCat,
  (dsId, xData, yData, zData, catData, xMeta, yMeta, zMeta, catMeta) => {
    if (!dsId || !xData || !yData || !zData || !catData) return undefined;
    let plotData = { id: "default", axes: {}, meta: {}, scale: {} };
    plotData.axes.x = xData || { values: [] };
    let xDomain = xMeta.xScale.domain().slice(0);
    let xmin = getQueryValue("xmin");
    if (xmin) {
      xDomain[0] = 1 * xmin;
    }
    let xmax = getQueryValue("xmax");
    if (xmax) {
      xDomain[1] = 1 * xmax;
    }
    if (xMeta.meta.clamp && xMeta.meta.clamp > xMeta.meta.limit[0]) {
      xDomain = clampDomain(xMeta.meta.limit[1], xMeta.meta.clamp, 25);
      xMeta.xScale.clamp(true);
    }
    xMeta.xScale.domain(xDomain);
    plotData.meta.x = xMeta;
    plotData.scale.x = xMeta.xScale;
    plotData.axes.y = yData || { values: [] };
    let yDomain = yMeta.xScale.domain().slice(0);
    let ymin = getQueryValue("ymin");
    if (ymin) {
      yDomain[0] = 1 * ymin;
    }
    let ymax = getQueryValue("ymax");
    if (ymax) {
      yDomain[1] = 1 * ymax;
    }
    if (yMeta.meta.clamp && yMeta.meta.clamp > yMeta.meta.limit[0]) {
      yDomain = clampDomain(yMeta.meta.limit[1], yMeta.meta.clamp, 25);
      // domain = [meta.clamp, range[1]]
      yMeta.xScale.clamp(true);
    }
    yMeta.xScale.domain(yDomain);
    plotData.meta.y = yMeta;
    plotData.scale.y = yMeta.xScale;
    plotData.axes.z = zData || { values: [] };
    plotData.meta.z = zMeta;
    plotData.scale.z = zMeta.xScale;
    plotData.axes.cat = catData;
    plotData.meta.cat = catMeta;
    plotData.scale.cat = catMeta.xScale;
    return plotData;
  }
);

export const getAllScatterPlotData = createSelector(
  getAllMainPlotData,
  getTransformFunction,
  (plotData, transform) => {
    if (!plotData) return undefined;
    let data = [];
    let scales = {};
    let axes = ["x", "y", "z"];
    if (
      plotData.axes.x.values.length == 0 ||
      plotData.axes.y.values.length == 0 ||
      plotData.axes.z.values.length == 0
    ) {
      return { data: [] };
    }
    axes.forEach((axis) => {
      scales[axis] = plotData.scale[axis].copy();
      if (plotData.scale[axis].hasOwnProperty("type")) {
        scales[axis].type = plotData.scale[axis].type;
      }
      // if (axis == 'z'){
      //   scales[axis] = d3.scaleSqrt().domain(scales[axis].domain())
      // }
      scales[axis].range([0, 900]);
    });
    let len = plotData.axes.x.values.length;
    let yClamp = plotData.meta.y.meta.clamp || Number.NEGATIVE_INFINITY;
    let yMin = plotData.meta.y.meta.limit[0];
    if (yClamp < yMin) {
      yClamp = Number.NEGATIVE_INFINITY;
    }
    let xClamp = plotData.meta.x.meta.clamp || Number.NEGATIVE_INFINITY;
    let xMin = plotData.meta.x.meta.limit[0];
    if (xClamp < xMin) {
      xClamp = Number.NEGATIVE_INFINITY;
    }
    for (let i = 0; i < len; i++) {
      let y = plotData.axes.y.values[i];
      let x = plotData.axes.x.values[i];
      if (
        x >= plotData.meta.x.meta.limit[0] &&
        x <= plotData.meta.x.meta.limit[1] &&
        y >= plotData.meta.y.meta.limit[0] &&
        y <= plotData.meta.y.meta.limit[1]
      ) {
        y = y < yClamp ? scales.y(yMin) : scales.y(y);
        x = x < xClamp ? scales.x(xMin) : scales.x(x);
        if (transform) [x, y] = transform([x, y]);
        data.push({
          id: i,
          x: x,
          y: 900 - y,
          z: plotData.axes.z.values[i],
        });
      }
    }
    return { data };
  }
);

export const getScatterPlotData = createSelector(
  getMainPlotData,
  getFilteredList,
  getTransformFunction,
  (plotData, list, transform) => {
    if (!plotData) return undefined;
    let data = [];
    let scales = {};
    let axes = ["x", "y", "z"];
    if (
      plotData.axes.x.values.length == 0 &&
      plotData.axes.y.values.length == 0 &&
      plotData.axes.z.values.length == 0
    ) {
      return { data: [] };
    }
    axes.forEach((axis) => {
      scales[axis] = plotData.meta[axis].xScale
        ? plotData.meta[axis].xScale.copy()
        : d3.scaleLinear().domain([0, 10]);
      if (axis == "z") {
        scales[axis] = d3.scaleSqrt().domain(scales[axis].domain());
        scales[axis].type = "scaleSqrt";
      } else {
        if (plotData.meta[axis].xScale.hasOwnProperty("type")) {
          scales[axis].type = plotData.meta[axis].xScale.type;
        }
      }
      scales[axis].range([0, 900]);
    });
    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;
    let len = Math.max(
      plotData.axes.x.values.length,
      plotData.axes.y.values.length,
      plotData.axes.z.values.length
    );
    let yClamp = Number.NEGATIVE_INFINITY;
    if (plotData.meta.y.meta) {
      let yClamp = plotData.meta.y.meta.clamp || Number.NEGATIVE_INFINITY;
      let yMin = plotData.meta.y.meta.limit[0];
      if (yClamp < yMin) {
        yClamp = Number.NEGATIVE_INFINITY;
      }
    }
    let xClamp = plotData.meta.x.meta.clamp || Number.NEGATIVE_INFINITY;
    let xMin = plotData.meta.x.meta.limit[0];
    if (xClamp < xMin) {
      xClamp = Number.NEGATIVE_INFINITY;
    }
    for (let i = 0; i < len; i++) {
      let y = plotData.axes.y.values[i];
      y = y < yClamp ? scales.y(yMin) : scales.y(y);
      let x = plotData.axes.x.values[i];
      x = x < xClamp ? scales.x(xMin) : scales.x(x);
      let z = plotData.axes.z.values[i];
      if (transform) [x, y] = transform([x, y]);
      data.push({
        id: list[i],
        index: i,
        x: x,
        y: 900 - y,
        z: z,
      });
      // if (i == 1){
      max = Math.max(max, z);
      min = Math.min(min, z);
      // }
    }
    let range = [min, max];
    return { data, range };
  }
);

export const getScatterPlotDataByCategory = createSelector(
  getMainPlotData,
  getScatterPlotData,
  getBinsForCat,
  getColorPalette,
  (plotData, scatterData, bins, palette) => {
    if (!plotData || !scatterData || !bins) return undefined;
    let data = [];
    let keys = {};
    if (plotData.axes.cat.values.length == 0) {
      return { data: [] };
    }
    bins.forEach((bin, i) => {
      data[i] = [];
      bin.keys.forEach((key) => {
        keys[key] = i;
      });
    });
    let len = plotData.axes.x.values.length;
    for (let i = 0; i < len; i++) {
      data[keys[plotData.axes.cat.values[i]]].push(scatterData.data[i]);
    }
    return { data, bins };
  }
);

const flattenNestedArray = (arr, errorBars, details) => {
  let func = (value, index) => value;
  // Experimental code to display relative values
  // if (reference && reference.values && details) {
  //   let scale = d3[details.meta.scale]()
  //     .domain(details.range)
  //     .clamp(true)
  //     .range([0, 1]);
  //   let midpoint = (scale(details.range[1]) + scale(details.range[0])) / 2;
  //   func = (value, index) => {
  //     return scale.invert(
  //       scale(value) - scale(reference.values[index]) + midpoint
  //     );
  //   };
  // }
  let newArr = [];
  let sd = [];
  let hasSd;
  for (let i = 0; i < arr.length; i++) {
    if (Array.isArray(arr[i])) {
      newArr[i] = [];
      sd[i] = [];
      for (let j = 0; j < arr[i].length; j++) {
        if (Array.isArray(arr[i][j])) {
          newArr[i][j] = func(arr[i][j][0], i);
          if (arr[i][j].length > 1) {
            if (errorBars) {
              hasSd = true;
              if (errorBars == "sd") {
                sd[i][j] = arr[i][j][1];
              } else if (errorBars == "se") {
                sd[i][j] = arr[i][j][1] / Math.sqrt(arr[i][j][2]);
              } else if (errorBars == "ci") {
                sd[i][j] = (1.96 * arr[i][j][1]) / Math.sqrt(arr[i][j][2]);
              }
            }
          }
        } else {
          newArr[i][j] = func(arr[i][j], i);
        }
      }
    } else {
      newArr[i] = func(arr[i], i);
    }
  }
  if (!hasSd) {
    sd = undefined;
  }
  return { values: newArr, sd };
};

const nestValues = (arr) => {
  let newArr = [];
  for (let i = 0; i < arr.length; i++) {
    newArr[i] = [arr[i]];
  }
  return newArr;
};

export const getWindowedXAxis = createSelector(
  getXAxis,
  getWindowSize,
  (axis, windowSize) => {
    let windowSuffix = windowSize == 0.1 ? "" : `_${windowSize}`;
    return `${axis}_windows${windowSuffix}`;
  }
);

export const getWindowedYAxis = createSelector(
  getYAxis,
  getWindowSize,
  (axis, windowSize) => {
    let windowSuffix = windowSize == 0.1 ? "" : `_${windowSize}`;
    return `${axis}_windows${windowSuffix}`;
  }
);

export const getWindowedZAxis = createSelector(
  getZAxis,
  getWindowSize,
  (axis, windowSize) => {
    let windowSuffix = windowSize == 0.1 ? "" : `_${windowSize}`;
    return `${axis}_windows${windowSuffix}`;
  }
);

export const getWindowedCatAxis = createSelector(
  getCatAxis,
  getWindowSize,
  (axis, windowSize) => {
    let windowSuffix = windowSize == 0.1 ? "" : `_${windowSize}`;
    return `${axis}_windows${windowSuffix}`;
  }
);

const getFilteredWindowDataForX = createSelector(
  getXAxis,
  getWindowedXAxis,
  getFields,
  (state) => getFilteredDataForFieldId(state, getWindowedXAxis(state)),
  getFilteredDataForX,
  getDetailsForX,
  getErrorBars,
  (axis, windowedAxis, fields, data, mainData, details, errorBars) => {
    if (fields[windowedAxis]) {
      if (data && data.values) {
        if (
          axis.match("proportion") ||
          axis.match("position") ||
          axis.match("length")
        ) {
          errorBars = false;
        }
        let { values, sd } = flattenNestedArray(
          data.values,
          errorBars,
          details
        );
        return { ...data, values, sd };
      }
      return data;
    }
    if (mainData && mainData.values) {
      let values = nestValues(mainData.values);
      return { ...mainData, values };
    }
    return mainData;
  }
);

const getFilteredWindowDataForY = createSelector(
  getYAxis,
  getWindowedYAxis,
  getFields,
  (state) => getFilteredDataForFieldId(state, getWindowedYAxis(state)),
  getFilteredDataForY,
  getDetailsForY,
  getErrorBars,
  (axis, windowedAxis, fields, data, mainData, details, errorBars) => {
    if (fields[windowedAxis]) {
      if (data && data.values) {
        if (
          axis.match("proportion") ||
          axis.match("position") ||
          axis.match("length")
        ) {
          errorBars = false;
        }
        let { values, sd } = flattenNestedArray(
          data.values,
          errorBars,
          details
        );
        return { ...data, values, sd };
      }
      return data;
    }
    if (mainData && mainData.values) {
      let values = nestValues(mainData.values);
      return { ...mainData, values };
    }
    return mainData;
  }
);

const getFilteredWindowDataForZ = createSelector(
  getZAxis,
  getWindowedZAxis,
  getFields,
  (state) => getFilteredDataForFieldId(state, getWindowedZAxis(state)),
  getFilteredDataForZ,
  getErrorBars,
  (axis, windowedAxis, fields, data, mainData, errorBars) => {
    if (fields[windowedAxis]) {
      if (data && data.values) {
        let { values, sd } = flattenNestedArray(data.values, errorBars);
        return { ...data, values, sd };
      }
      return data;
    }
    if (mainData && mainData.values) {
      let values = nestValues(mainData.values);
      return { ...mainData, values };
    }
    return mainData;
  }
);

const getFilteredWindowDataForCat = createSelector(
  getCatAxis,
  getWindowedCatAxis,
  getFields,
  (state) => getFilteredDataForFieldId(state, getWindowedCatAxis(state)),
  getFilteredDataForCat,
  (axis, windowedAxis, fields, data, mainData) => {
    if (fields[windowedAxis]) {
      if (data && data.values) {
        let { values } = flattenNestedArray(data.values);
        return { ...data, values };
      }
      return data;
    }
    if (mainData && mainData.values) {
      let values = nestValues(mainData.values);
      return { ...mainData, values };
    }
    return mainData;
  }
);

export const getWindowPlotData = createSelector(
  getDatasetID,
  getFilteredWindowDataForX,
  getFilteredWindowDataForY,
  getFilteredWindowDataForZ,
  getFilteredWindowDataForCat,
  getFilteredDataForCat,
  getDetailsForX,
  getDetailsForY,
  getDetailsForZ,
  getDetailsForCat,
  (
    dsId,
    xData,
    yData,
    zData,
    catData,
    mainCatData,
    xMeta,
    yMeta,
    zMeta,
    catMeta
  ) => {
    if (!dsId || !catData) return undefined;
    let plotData = { id: "default", axes: {}, meta: {} };
    plotData.axes.x = xData || { values: [] };
    let xDomain = xMeta.xScale ? xMeta.xScale.domain().slice(0) : [0, 10];
    let xmin = getQueryValue("xmin");
    if (xmin) {
      xDomain[0] = 1 * xmin;
    }
    let xmax = getQueryValue("xmax");
    if (xmax) {
      xDomain[1] = 1 * xmax;
    }
    xMeta.xScale.domain(xDomain);
    xMeta.xScale.type = xMeta.meta.scale;
    plotData.meta.x = xMeta;
    plotData.axes.y = yData || { values: [] };
    let yDomain = yMeta.xScale ? yMeta.xScale.domain().slice(0) : [0, 10];
    let ymin = getQueryValue("ymin");
    if (ymin) {
      yDomain[0] = 1 * ymin;
    }
    let ymax = getQueryValue("ymax");
    if (ymax) {
      yDomain[1] = 1 * ymax;
    }
    if (yMeta.xScale) {
      yMeta.xScale.domain(yDomain);
      yMeta.xScale.type = yMeta.meta.scale;
    }
    plotData.meta.y = yMeta;
    plotData.axes.z = zData || { values: [] };
    plotData.meta.z = zMeta;
    plotData.axes.cat = catData;
    plotData.axes.mainCat = mainCatData;
    plotData.meta.cat = catMeta;
    return plotData;
  }
);

export const getWindowBinsForCat = createSelector(
  getBinsForCat,
  getFilteredWindowDataForCat,
  getOtherLimit,
  getColorPalette,
  (bins, cats, otherLimit, palette) => {
    if (!bins || !cats) {
      return bins;
    }
    let keys = {};
    let windowBins = [];
    let others = [];
    bins.forEach((bin, i) => {
      bin.keys.forEach((key) => {
        if (!keys[bin.id]) {
          keys[bin.id] = [];
        }
        keys[bin.id].push(key);
      });
      windowBins.push({ ...bin });
    });
    cats.keys.forEach((id, key) => {
      if (!keys[id]) {
        if (windowBins.length < otherLimit) {
          if (
            windowBins.length == otherLimit - 1 &&
            cats.keys.length > otherLimit
          ) {
            id = "other";
          }
          windowBins.push({
            id,
            keys: [key],
            x1: key,
            x2: key + 1,
            color: palette.colors[windowBins.length - 1],
          });
        } else {
          others.push(key);
        }
      }
    });
    if (others.length > 0) {
      windowBins[windowBins.length - 1].keys.push(...others);
    }
    return windowBins;
  }
);

export const getLinesPlotData = createSelector(
  getWindowPlotData,
  getFilteredList,
  getTransformFunction,
  getDetailsForZ,
  getZScale,
  getPlotResolution,
  getPlotScale,
  getWindowBinsForCat,
  getColorPalette,
  (
    plotData,
    list,
    transform,
    details,
    scale,
    res,
    plotScale,
    bins,
    palette
  ) => {
    if (!plotData) return {};
    let coords = [];
    let scales = {};
    let axes = ["x", "y", "z"];
    if (
      plotData.axes.x.values.length == 0 ||
      plotData.axes.y.values.length == 0 ||
      plotData.axes.z.values.length == 0 ||
      plotData.axes.cat.values.length == 0
    ) {
      return { coords: [] };
    }
    let zScale = d3[scale]()
      .domain(details.range)
      .range([2, (2 * 900) / res]);
    axes.forEach((axis) => {
      scales[axis] = plotData.meta[axis].xScale
        ? plotData.meta[axis].xScale.copy()
        : d3.scaleLinear().domain([0, 10]);
      if (axis == "z") {
        scales[axis] = d3.scaleSqrt().domain(scales[axis].domain());
        scales[axis].type = "scaleSqrt";
      } else {
        if (plotData.meta[axis].xScale.hasOwnProperty("type")) {
          scales[axis].type = plotData.meta[axis].xScale.type;
        }
      }
      scales[axis].range([0, 900]);
    });
    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;
    let len = Math.max(
      plotData.axes.x.values.length,
      plotData.axes.y.values.length,
      plotData.axes.z.values.length
    );
    let yClamp = Number.NEGATIVE_INFINITY;
    if (plotData.meta.y.meta) {
      let yClamp = plotData.meta.y.meta.clamp || Number.NEGATIVE_INFINITY;
      let yMin = plotData.meta.y.meta.limit[0];
      if (yClamp < yMin) {
        yClamp = Number.NEGATIVE_INFINITY;
      }
    }
    let xClamp = plotData.meta.x.meta.clamp || Number.NEGATIVE_INFINITY;
    let xMin = plotData.meta.x.meta.limit[0];
    if (xClamp < xMin) {
      xClamp = Number.NEGATIVE_INFINITY;
    }
    let keys = {};
    bins.forEach((bin, i) => {
      bin.keys.forEach((key) => {
        keys[key] = i;
      });
    });
    axes = ["x", "y", "z", "cat"];
    for (let i = 0; i < len; i++) {
      let xs = [];
      let ys = [];
      let zs = [];
      let rs = [];
      let cats = [];
      let xsd = [];
      let ysd = [];
      let wins = Math.max(
        plotData.axes.x.values[i].length,
        plotData.axes.y.values[i].length,
        plotData.axes.z.values[i].length
        // plotData.axes.cat.values[i].length
      );
      let zSize = -1;
      for (let j = 0; j < wins; j++) {
        let values = {};
        let sd = {};
        let valid = true;
        axes.forEach((axis) => {
          if (plotData.axes[axis].values[i].length == wins) {
            values[axis] = plotData.axes[axis].values[i][j];
            if (plotData.axes[axis].sd) {
              sd[axis] = plotData.axes[axis].sd[i][j];
            }
          } else {
            values[axis] = plotData.axes[axis].values[i][0];
            if (plotData.axes[axis].sd) {
              sd[axis] = plotData.axes[axis].sd[i][0];
            }
          }
          if (isNaN(values[axis]) || values[axis] === undefined) {
            valid = false;
          }
        });
        if (!valid) {
          continue;
        }
        // ensure last bin is not too small to plot
        // if (values.z >= zSize) {
        //   zSize = values.z;
        // } else {
        //   if (values.z < zSize * 0.75) {
        //     break;
        //   }
        // }
        values.y = values.y < yClamp ? scales.y(yMin) : scales.y(values.y);
        values.x = values.x < xClamp ? scales.x(xMin) : scales.x(values.x);
        values.cat = keys[values.cat];
        if (transform) [values.x, values.y] = transform([values.x, values.y]);
        xs.push(values.x);
        ys.push(900 - values.y);
        if (sd.x) {
          sd.x = scales.x(scales.x.domain()[0] + sd.x);
          xsd.push(sd.x);
        }
        if (sd.y) {
          sd.y =
            scales.y(scales.y.domain()[0] + sd.y) -
            scales.y(scales.y.domain()[0]);
          ysd.push(sd.y);
        }
        zs.push(values.z);
        rs.push(zScale(values.z) * plotScale);
        cats.push(values.cat);
        max = Math.max(max, values.z);
        min = Math.min(min, values.z);
      }
      if (xsd.length == 0) {
        xsd = undefined;
      }
      if (ysd.length == 0) {
        ysd = undefined;
      }
      coords.push({
        id: list[i],
        index: i,
        x: xs,
        y: ys,
        z: zs,
        r: rs,
        cats: cats,
        cat: keys[plotData.axes.mainCat.values[i]],
        xsd,
        ysd,
      });
    }
    let range = [min, max];
    return { coords, range, colors: palette.colors };
  }
);
// const weightedMeanSd = (values, weights, sumWeight, log) => {
//   if (log){
//     let sum = 0
//     let len = values.length
//     for (let i = 0; i < len; i++){
//       sum += Math.log10(values[i])*weights[i]
//     }
//     let mean = sum/sumWeight
//     let variance = 0
//     for (let i = 0; i < len; i++){
//       variance += weights[i] * ((Math.log10(values[i]) - mean) ** 2)
//     }
//     variance /= sumWeight
//     let stDev = Math.sqrt(variance)
//     let upper = 10 ** (mean + 2 * stDev)
//     let lower = 10 ** (mean - 2 * stDev)
//     mean = 10 ** mean
//     return {mean, stDev, upper, lower}
//   }
//   else {
//     let sum = 0
//     let len = values.length
//     for (let i = 0; i < len; i++){
//       sum += values[i]*weights[i]
//     }
//     let mean = sum/sumWeight
//     let variance = 0
//     for (let i = 0; i < len; i++){
//       variance += weights[i] * (values[i] - mean) ** 2
//     }
//     variance /= sumWeight
//     let stDev = Math.sqrt(variance)
//     let upper = mean + 2 * stDev
//     let lower = mean - 2 * stDev
//     return {mean, stDev, upper, lower}
//   }
// }

const weightedMedian = (arr, total) => {
  if (arr.length == 0) return undefined;
  let sum = 0;
  let mid = total / 2;
  let median;
  let sorted = arr.sort((a, b) => a.value - b.value);
  for (let i = 0; i < arr.length; i++) {
    sum += sorted[i].weight;
    if (sum > mid) {
      return sorted[i].value;
    }
  }
};

const calculateSlope = (xValues, yValues, zValues, xLog, yLog) => {
  let sum = 0,
    xSum = 0,
    ySum = 0,
    sumSq = 0;
  let n = xValues.length;
  for (let i = 0; i < n; i++) {
    let x = xValues[i];
    let y = yValues[i];
    sum += x * y;
    xSum += x;
    ySum += y;
    sumSq += x * x;
  }
  let a = sum * n;
  let b = xSum * ySum;
  let c = sumSq * n;
  let d = xSum * xSum;
  let m = (a - b) / (c - d);
  return m;
};

function weightedLinearRegression(xValues, yValues, weights) {
  let sums = { w: 0, wx: 0, wx2: 0, wy: 0, wxy: 0 };
  let n = xValues.length;
  for (let i = 0; i < n; i++) {
    let x = xValues[i];
    let y = yValues[i];
    let w = weights[i];
    sums.w += w;
    sums.wx += x * w;
    sums.wx2 += x * x * w;
    sums.wy += y * w;
    sums.wxy += x * y * w;
  }
  let denominator = sums.w * sums.wx2 - sums.wx * sums.wx;
  let m = (sums.w * sums.wxy - sums.wx * sums.wy) / denominator;
  let c = (sums.wy * sums.wx2 - sums.wx * sums.wxy) / denominator;
  return { m, c };
}

const weightedMeanSd = (values, weights, sumWeight, log) => {
  log = false;
  if (log) {
    let sum = 0;
    let len = values.length;
    let arr = [];
    for (let i = 0; i < len; i++) {
      sum += Math.log10(values[i]) * weights[i];
      arr.push({ weight: weights[i], i, value: values[i] });
    }
    let median = weightedMedian(arr, sumWeight);
    let mean = sum / sumWeight;
    let variance = 0;
    for (let i = 0; i < len; i++) {
      variance += weights[i] * (Math.log10(values[i]) - mean) ** 2;
    }
    variance /= sumWeight;
    let stDev = Math.sqrt(variance);
    let upper = 10 ** (mean + 2 * stDev);
    let lower = 10 ** (mean - 2 * stDev);
    mean = 10 ** mean;
    return { mean, median, stDev, upper, lower };
  } else {
    let sum = 0;
    let len = values.length;
    let arr = [];
    for (let i = 0; i < len; i++) {
      sum += values[i] * weights[i];
      arr.push({ weight: weights[i], i, value: values[i] });
    }
    let median = weightedMedian(arr, sumWeight);
    let mean = sum / sumWeight;
    let variance = 0;
    for (let i = 0; i < len; i++) {
      variance += weights[i] * (values[i] - mean) ** 2;
    }
    variance /= sumWeight;
    let stDev = Math.sqrt(variance);
    let upper = mean + 2 * stDev;
    let lower = mean - 2 * stDev;
    return { mean, median, stDev, upper, lower };
  }
};

export const scaleFactor = createSelector(getMainPlotData, (plotData) => {
  let yScale = plotData.meta.y.xScale.copy().clamp(false);
  let xScale = plotData.meta.x.xScale.copy().clamp(false);
  yScale.range([0, 900]);
  xScale.range([0, 900]);
  if (plotData.meta.x.xScale.hasOwnProperty("type")) {
    xScale.type = plotData.meta.x.xScale.type;
  }
  if (plotData.meta.y.xScale.hasOwnProperty("type")) {
    yScale.type = plotData.meta.y.xScale.type;
  }
  let yRange = 900;
  let yDenom = yScale.domain()[1] - yScale.domain()[0];
  if (yScale.type == "scaleLog") {
    yDenom = Math.log10(yDenom);
  } else if (yScale.type == "scaleSqrt") {
    yDenom = Math.sqrt(yDenom);
  }
  let yFactor = yRange / yDenom;
  let xRange = 900;
  let xDenom = xScale.domain()[1] - xScale.domain()[0];
  if (xScale.type == "scaleLog") {
    xDenom = Math.log10(xDenom);
  } else if (xScale.type == "scaleSqrt") {
    xDenom = Math.sqrt(xDenom);
  }
  let xFactor = xRange / xDenom;
  let factor = xFactor / yFactor;
  return { factor, xScale, yScale };
});

export const getKitePlotData = createSelector(
  getMainPlotData,
  getBinsForCat,
  getTransformFunction,
  getColorPalette,
  (plotData, bins, transform, palette) => {
    if (!plotData || !bins) return undefined;
    let data = [];
    let keys = {};
    if (plotData.axes.cat.values.length == 0) {
      return { data: [] };
    }
    bins.forEach((bin, i) => {
      data[i] = { x: [], y: [], z: [], zSum: 0 };
      bin.keys.forEach((key) => {
        keys[key] = i;
      });
    });
    let yScale = plotData.meta.y.xScale.copy().clamp(false);
    let xScale = plotData.meta.x.xScale.copy().clamp(false);
    if (plotData.meta.y.xScale.hasOwnProperty("type")) {
      yScale.type = plotData.meta.y.xScale.type;
    }
    if (plotData.meta.x.xScale.hasOwnProperty("type")) {
      xScale.type = plotData.meta.x.xScale.type;
    }
    yScale.range([0, 900]);
    xScale.range([0, 900]);
    let len = plotData.axes.x.values.length;
    for (let i = 0; i < len; i++) {
      let bin = keys[plotData.axes.cat.values[i]];
      let x = plotData.axes.x.values[i];
      let y = plotData.axes.y.values[i];
      let z = plotData.axes.z.values[i];
      if (!plotData.meta.x.meta.clamp || x >= plotData.meta.x.meta.clamp) {
        if (!plotData.meta.y.meta.clamp || y >= plotData.meta.y.meta.clamp) {
          x = xScale(x);
          y = yScale(y);
          if (transform) [x, y] = transform([x, y]);
          data[bin].x.push(x);
          data[bin].y.push(900 - y);
          data[bin].z.push(z);
          data[bin].zSum += z;
        }
      }
    }
    let coords = [],
      equations = [];
    data.forEach((bin, i) => {
      coords[i] = {};
      equations[i] = { order: 1, x: 0 };
      bin.xWeighted = weightedMeanSd(
        bin.x,
        bin.z,
        bin.zSum,
        plotData.meta.x.meta.scale == "scaleLog"
      );
      bin.yWeighted = weightedMeanSd(
        bin.y,
        bin.z,
        bin.zSum,
        plotData.meta.y.meta.scale == "scaleLog"
      );
      if (bin.xWeighted.lower) {
        coords[i].x = [
          [bin.xWeighted.lower, bin.yWeighted.median],
          [bin.xWeighted.upper, bin.yWeighted.median],
        ];
        coords[i].y = [
          [bin.xWeighted.median, bin.yWeighted.lower],
          [bin.xWeighted.median, bin.yWeighted.upper],
        ];
        coords[i].poly = [
          [bin.xWeighted.lower, bin.yWeighted.median],
          [bin.xWeighted.median, bin.yWeighted.lower],
          [bin.xWeighted.upper, bin.yWeighted.median],
          [bin.xWeighted.median, bin.yWeighted.upper],
          [bin.xWeighted.lower, bin.yWeighted.median],
        ];
        coords[i].weight = Math.ceil(Math.log10(bin.x.length));
        if (bin.x.length >= 30) {
          let eqn = weightedLinearRegression(bin.x, bin.y, bin.z);
          let angle = (Math.tan(eqn.m) * 180) / Math.PI;
          coords[i].angle = angle;
          coords[i].radians = eqn.m;
          equations[i].factor = eqn.m;
          equations[i].intercept = -eqn.m * bin.xWeighted.median;
          equations[i].origin = {
            x: bin.xWeighted.median,
            y: bin.yWeighted.median,
          };
          equations[i].index = i;
        } else {
          coords[i].angle = 0;
          coords[i].radians = 0;
        }
      }
    });

    return { coords, equations, bins, xScale, yScale, colors: palette.colors };
  }
);

export const getTransformEquation = createSelector(
  getTransformFunctionParams,
  scaleFactor,
  (params, scaleFactor) => {
    let factor = params.factor * scaleFactor.factor;
    let intercept;
    if (params.intercept == 0) {
      intercept = 0;
    } else if (scaleFactor.yScale.type == "scaleLog") {
      intercept = Math.log10(scaleFactor.yScale.invert(params.intercept));
    } else {
      intercept = scaleFactor.yScale.invert(params.intercept);
    }
    let xType = scaleFactor.xScale.type;
    let yType = scaleFactor.yScale.type;
    return { factor, intercept, xType, yType };
  }
);
// export const getKitePlotData = createSelector(
//   getMainPlotData,
//   getBinsForCat,
//   getColorPalette,
//   (plotData,bins,palette) => {
//     if (!plotData || !bins) return undefined
//     let data = [];
//     let keys = {}
//     if (plotData.axes.cat.values.length == 0){
//       return {data:[]}
//     }
//     bins.forEach((bin,i)=>{
//       data[i] = {x:[],y:[],z:[],zSum:0}
//       bin.keys.forEach(key=>{
//         keys[key] = i
//       })
//     })
//     let len = plotData.axes.x.values.length
//     for (let i = 0; i < len; i++){
//       let bin = keys[plotData.axes.cat.values[i]]
//       let x = plotData.axes.x.values[i]
//       let y = plotData.axes.y.values[i]
//       let z = plotData.axes.z.values[i]
//       if (!plotData.meta.x.meta.clamp || x >= plotData.meta.x.meta.clamp){
//         if (!plotData.meta.y.meta.clamp || y >= plotData.meta.y.meta.clamp){
//           data[bin].x.push(x)
//           data[bin].y.push(y)
//           data[bin].z.push(z)
//           data[bin].zSum += z
//         }
//       }
//     }
//     let coords = []
//     let yScale = plotData.meta.y.xScale.copy()
//     let xScale = plotData.meta.x.xScale.copy()
//     yScale.range([900,0])
//     xScale.range([0,900])
//     data.forEach((bin,i)=>{
//       coords[i] = {}
//       bin.xWeighted = weightedMeanSd(bin.x, bin.z, bin.zSum, plotData.meta.x.meta.scale == 'scaleLog')
//       bin.yWeighted = weightedMeanSd(bin.y, bin.z, bin.zSum, plotData.meta.y.meta.scale == 'scaleLog')
//       if (bin.xWeighted.lower){
//         coords[i].x = [
//           [xScale(bin.xWeighted.lower),yScale(bin.yWeighted.median)],
//           [xScale(bin.xWeighted.upper),yScale(bin.yWeighted.median)]
//         ]
//         coords[i].y = [
//           [xScale(bin.xWeighted.median),yScale(bin.yWeighted.lower)],
//           [xScale(bin.xWeighted.median),yScale(bin.yWeighted.upper)]
//         ]
//         coords[i].poly = [
//           [xScale(bin.xWeighted.lower),yScale(bin.yWeighted.median)],
//           [xScale(bin.xWeighted.median),yScale(bin.yWeighted.lower)],
//           [xScale(bin.xWeighted.upper),yScale(bin.yWeighted.median)],
//           [xScale(bin.xWeighted.median),yScale(bin.yWeighted.upper)],
//           [xScale(bin.xWeighted.lower),yScale(bin.yWeighted.median)]
//         ]
//         if (bin.x.length >= 30){
//           let m = calculateSlope(bin.x, bin.y, plotData.meta.x.meta.scale == 'scaleLog', plotData.meta.y.meta.scale == 'scaleLog')
//
//           coords[i].slope = m
//         }
//         else {
//           coords[i].slope = 0
//         }
//       }
//     })
//     return {coords,bins,colors:palette.colors};
//   }
// )

const sliceObject = (obj, index) => {
  let slice = {};
  Object.keys(obj).forEach((key) => {
    slice[key] = obj[key].slice(index, index + 1)[0];
  });
  return slice;
};
const createSelectorForCategoryIndex = byIdSelectorCreator();
const createSelectorForCircleCategoryIndex = byIdSelectorCreator();
const createSelectorForSquareCategoryIndex = byIdSelectorCreator();
const createSelectorForSliceIndex = byIdSelectorCreator();
const createSelectorForColorIndex = byIdSelectorCreator();

const _getCategoryIndexAsMemoKey = (state, categoryIndex) => categoryIndex;
//const getScatterPlotDataForCategory = (state, categoryIndex) => getScatterPlotDataByCategory(state)[categoryIndex];
const getColorByIndex = createSelectorForColorIndex(
  _getCategoryIndexAsMemoKey,
  _getCategoryIndexAsMemoKey,
  getColorPalette,
  (index, palette) => {
    return palette.colors[index];
  }
);

export const getScatterPlotDataSlice = createSelectorForSliceIndex(
  _getCategoryIndexAsMemoKey,
  _getCategoryIndexAsMemoKey,
  getScatterPlotDataByCategory,
  (index, plotData) => {
    plotData = sliceObject(plotData, index);
    return plotData;
  }
);

export const getScatterPlotDataForCategoryIndex =
  createSelectorForCategoryIndex(
    _getCategoryIndexAsMemoKey,
    getScatterPlotDataSlice,
    getColorByIndex,
    (plotData, color, res) => {
      plotData.color = color;
      return plotData;
    }
  );

export const getCirclePlotDataForCategoryIndex =
  createSelectorForCircleCategoryIndex(
    _getCategoryIndexAsMemoKey,
    getScatterPlotDataForCategoryIndex,
    getDetailsForZ,
    getZScale,
    getPlotResolution,
    getPlotScale,
    (catData, details, scale, res, plotScale) => {
      let zScale = d3[scale]()
        .domain(details.range)
        .range([2, (2 * 900) / res]);
      if (catData.data) {
        catData.data.forEach((datum) => {
          datum.r = zScale(datum.z) * plotScale;
        });
      } else {
        catData.data = [];
      }
      return catData;
    }
  );

export const getSquareBinPlotDataForCategoryIndex =
  createSelectorForSquareCategoryIndex(
    _getCategoryIndexAsMemoKey,
    getScatterPlotDataForCategoryIndex,
    (plotData) => {
      let size = 900; // FIXME: magic number
      let res = 20; // FIXME: magic number
      let side = size / res;
      let squares = [];
      for (let i = 0; i <= res; i++) {
        squares[i] = [];
        for (let j = 0; j <= res; j++) {
          squares[i][j] = {
            id: i + "_" + j,
            x: i,
            y: j,
            ids: [],
            zs: [],
            indices: [],
          };
        }
      }
      plotData.data.forEach((datum) => {
        let x = Math.floor(datum.x / side);
        let y = Math.floor(datum.y / side);
        squares[x][y].ids.push(datum.id);
        squares[x][y].indices.push(datum.index);
        squares[x][y].zs.push(datum.z);
      });
      let data = [];
      for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
          if (squares[i][j].ids.length > 0) {
            data.push(squares[i][j]);
          }
        }
      }
      plotData.data = data;
      plotData.side = side;
      return plotData;
    }
  );

const createSelectorForMainPlotCategory = byIdSelectorCreator();

export const getCategoryListForMainPlot = createSelectorForMainPlotCategory(
  getAxisTitle,
  (state) => getPlainCategoryListForFieldId(state, getAxisTitle(state, "cat")),
  (list) => {
    return list;
  }
);
