import * as d3 from "d3";

import { createAction, handleAction, handleActions } from "redux-actions";
import {
  getBinsForFieldId,
  getDetailsForFieldId,
  getRawDataForFieldId,
  histogram,
} from "./field";
import { getDimensionsbyDimensionId, getPreviewDimensions } from "./dimension";
import { getFilteredList, updateFilterList } from "./filter";

import { byIdSelectorCreator } from "./selectorCreators";
import { createSelector } from "reselect";
import deep from "deep-get-set";
import { getColorPalette } from "./color";
import { getMainPlot } from "./plot";
import { getOtherLimit } from "./plotParameters";
import { getQueryValue } from "./location";
import { getSelectedDatasetMeta } from "./dataset";
import immutableUpdate from "immutable-update";
import store from "../store";

const _getFieldIdAsMemoKey = (state, fieldId) => fieldId;
const _getFilterIdAsMemoKey = (state, filterId) => filterId;

const createSelectorForFilterId = byIdSelectorCreator();
const getMetaDataForFilter = (state, filterId) =>
  state.filters ? state.filters.byId[filterId] : {};

export const getDetailsForFilterId = createSelectorForFilterId(
  _getFilterIdAsMemoKey,
  getMetaDataForFilter,
  getDetailsForFieldId,
  getBinsForFieldId,
  (filterMeta = {}, fieldMeta = {}, bins) => {
    let obj = {
      filterId: filterMeta.id,
    };
    if (fieldMeta.meta.type == "variable") {
      obj.filterType = "range";
      let range = fieldMeta.range ? fieldMeta.range.slice(0) : [1, 10];
      range = [
        1 * getQueryValue("min" + filterMeta.id) || range[0],
        1 * getQueryValue("max" + filterMeta.id) || range[1],
      ];
      if (isNaN(range[0])) range[0] = 0;
      if (isNaN(range[1])) range[1] = 10;
      obj.filterRange = filterMeta.range ? filterMeta.range.slice(0) : range;
      if (isNaN(obj.filterRange[0])) obj.filterRange[0] = 0;
      if (isNaN(obj.filterRange[1])) obj.filterRange[1] = 10;
      obj.filterLimit = range;
      obj.fieldLimit = fieldMeta.meta.limit || range;
      obj.xScale = fieldMeta.xScale;
    }
    if (fieldMeta.meta.type == "category") {
      obj.filterType = "category";
      obj.keys = filterMeta.keys;
      obj.toggled = [];
      if (obj.keys.length > 0) {
        bins.forEach((bin, i) => {
          obj.toggled[i] = bin.keys.some((key) => {
            return obj.keys.indexOf(key) != -1;
          });
        });
      } else {
        //obj.toggled = filterMeta.toggled
      }
    }
    if (fieldMeta.meta.type == "selection") {
      obj.filterType = "selection";
    }
    obj.clonedFrom = fieldMeta.clonedFrom;
    obj.invert = filterMeta.invert;
    return obj;
  }
);

export const getRawDataForLength = createSelector(
  (state) => getRawDataForFieldId(state, "length"),
  (data) => data
);

export const getRawDataForNCount = createSelector(
  (state) => getRawDataForFieldId(state, "ncount"),
  (data) => data
);

// const createAdjustedDataSelectorForFieldId = byIdSelectorCreator();

// export const getAdjustedDataForFieldId = createAdjustedDataSelectorForFieldId(
//   _getFieldIdAsMemoKey,
//   _getFieldIdAsMemoKey,
//   getRawDataForFieldId,
//   getRawDataForLength,
//   getRawDataForNCount,
//   getAdjustCoverage,
//   (fieldId, rawData, length, ncount, adjustCoverage) => {
//     if (!rawData) return undefined
//     if (!adjustCoverage) return rawData
//     if (!length || !ncount) return rawData
//     if (fieldId.match(/_cov$/)){
//       let values = []
//       if (rawData.values){
//         let raw = rawData.values
//         let lengths = length.values
//         let counts = ncount.values
//         let len = raw.length
//         for (var i = 0; i < len; i++){
//           values.push(raw[i]/(1-(counts[i]/lengths[i])));
//         }
//       }
//       return {values,keys:rawData.keys}
//     }
//     else {
//       return rawData
//     }
//   }
// );

const createFilteredDataSelectorForFieldId = byIdSelectorCreator();

export const getFilteredDataForFieldId = createFilteredDataSelectorForFieldId(
  _getFieldIdAsMemoKey,
  getFilteredList,
  getRawDataForFieldId,
  (list = [], rawData) => {
    if (!rawData) return undefined;
    let values = [];
    if (rawData.values) {
      let raw = rawData.values;
      let len = list.length;
      for (var i = 0; i < len; i++) {
        values.push(raw[list[i]]);
      }
    }
    return { values, keys: rawData.keys };
  }
);

const createFilteredBarSelectorForFieldId = byIdSelectorCreator();

export const getFilteredBarsForFieldId = createFilteredBarSelectorForFieldId(
  _getFieldIdAsMemoKey,
  getFilteredDataForFieldId,
  getBinsForFieldId,
  getDetailsForFieldId,
  getPreviewDimensions,
  getOtherLimit,
  (data, fieldBins = [], details = {}, dimensions, otherLimit) => {
    let bars = [];
    if (data) {
      let x = details.xScale;
      let bins = [];
      if (details.meta.type == "variable") {
        x.range([0, 25]);
        let thresh = Array.from(Array(24).keys()).map((n) => {
          return x.invert(n + 1);
        });
        if (details.meta.clamp && details.meta.clamp > details.meta.limit[0]) {
          bins = histogram(x, thresh.slice(0), data.values, details.meta.clamp);
        } else {
          bins = d3.histogram().domain(x.domain()).thresholds(thresh)(
            data.values
          );
        }
      }
      if (details.meta.type == "category") {
        let nested = d3
          .nest()
          .key((d) => data.keys[d])
          .rollup((d) => d.length)
          .entries(data.values);
        fieldBins.forEach((obj, i) => {
          let index = nested.findIndex((n) => obj.id == n.key);
          if (index > -1) {
            bins[i] = {
              id: obj.id,
              x0: i,
              x1: i + 1,
              length: nested[index].value,
            };
          } else {
            bins[i] = { id: obj.id, x0: i, x1: i + 1, length: 0 };
          }
        });
        if (nested.length > otherLimit) {
          let sortedSum = bins.map((a) => a.length).reduce((a, b) => a + b);
          let nestedSum = nested.map((a) => a.value).reduce((a, b) => a + b);
          let other = nestedSum - sortedSum;
          bins[otherLimit - 1].length = nestedSum - sortedSum;
        }
      }

      let y = d3
        .scaleLinear()
        .domain([
          0,
          d3.max(fieldBins, function (d) {
            return d.length;
          }),
        ])
        .range([dimensions.height, 0]);
      if (details.meta.type == "category") {
        x = d3.scaleLinear();
        x.domain([0, 10]);
        y = d3
          .scaleSqrt()
          .domain([
            0,
            d3.max(fieldBins, function (d) {
              return d.length;
            }),
          ])
          .range([dimensions.height, 0]);
      }
      x.range([0, dimensions.width]);
      bins.forEach((d, i) => {
        bars.push({
          id: d.id || i,
          x: x(d.x0),
          y: y(d.length),
          width: x(d.x1) - x(d.x0) - 1,
          height: dimensions.height - y(d.length) || 0,
        });
      });
    }
    return bars;
  }
);

const createPlainCategoryListSelectorForFieldId = byIdSelectorCreator();
export const getPlainCategoryListForFieldId = createPlainCategoryListSelectorForFieldId(
  _getFieldIdAsMemoKey,
  _getFieldIdAsMemoKey,
  getBinsForFieldId,
  getColorPalette,
  getMainPlot,
  (fieldId, bins = [], palette, plot) => {
    bins.forEach((b, i) => {
      b.color =
        fieldId == plot.axes.cat ? palette.colors[i] : "rgb(215, 205, 204)";
    });
    return { bins, plot };
  }
);

const createCategoryListSelectorForFieldId = byIdSelectorCreator();

export const getCategoryListForFieldId = createCategoryListSelectorForFieldId(
  _getFieldIdAsMemoKey,
  getBinsForFieldId,
  getDetailsForFilterId, // FIXME
  getColorPalette,
  getMainPlot,
  (bins = [], filter = {}, palette, plot) => {
    let matchId = filter.filterId;
    // let matchId = filter.clonedFrom || filter.filterId
    bins.forEach((b, i) => {
      b.toggled = filter.toggled[i];
      b.color =
        plot.axes.cat == matchId ? palette.colors[i] : "rgb(215, 205, 204)";
      // b.color =  plot.axes.cat.match(filter.filterID) ? palette.colors[i] : 'rgb(215, 205, 204)'
    });
    return { bins, filter, plot };
  }
);

export const getFilteredSummary = createSelector(
  getSelectedDatasetMeta,
  getFilteredList,
  (meta = {}, list = []) => {
    let toPercent = d3.format(",.1%");
    return {
      count: meta.records || 0,
      type: meta.record_type || "",
      selected: list.length || 0,
      percentage: toPercent(list.length / meta.records || 0),
    };
  }
);
