import * as d3 from "d3";

import {
  getAllScatterPlotData,
  getFilteredDataForCat,
  getScatterPlotData,
} from "./plotData";
import { getBinsForCat, getBinsForFieldId } from "./field";
import { getDetailsForFieldId, getRawDataForFieldId } from "./field";
import {
  getPlotResolution,
  getPlotScale,
  getSideMax,
  getZReducer,
  getZScale,
} from "./plotParameters";

import { byIdSelectorCreator } from "./selectorCreators";
import { createSelector } from "reselect";
import { getColorPalette } from "./color";
import { getFilteredDataForFieldId } from "./preview";
import { getMainPlot } from "./plot";
import { getSelectedRecordsAsObject } from "./select";
import immutableUpdate from "immutable-update";

export const getSquareGrid = createSelector(getPlotResolution, (res) => {
  let size = 900.1; // FIXME: magic number
  let height = size / (res * 1 + 1);
  let width = height;
  let data = [];
  for (let i = 0; i <= res; i++) {
    for (let j = 0; j <= res; j++) {
      data.push({
        id: i + "_" + j,
        i,
        j,
        x: i * width,
        y: j * height,
        height,
        width,
      });
    }
  }
  return { data, size, res, height, width };
});

export const setCoords = (datum, grid) => {
  let coords = [
    Math.floor(datum.x / grid.width),
    Math.floor(datum.y / grid.height),
  ];
  return coords;
};

export const getAllScatterPlotDataBySquareBin = createSelector(
  getSquareGrid,
  getAllScatterPlotData,
  (grid, scatterData) => {
    if (!scatterData) return undefined;
    let squares = [];
    grid.data.forEach((d) => {
      squares[d.i] = squares[d.i] || [];
      squares[d.i][d.j] = {
        id: d.id,
        i: d.i,
        j: d.j,
        ids: [],
        zs: [],
        indices: [],
      };
    });
    scatterData.data.forEach((datum) => {
      let coords = setCoords(datum, grid);
      if (squares[coords[0]] && squares[coords[0]][coords[1]]) {
        squares[coords[0]][coords[1]].ids.push(datum.id);
        squares[coords[0]][coords[1]].indices.push(datum.index);
        squares[coords[0]][coords[1]].zs.push(datum.z);
      }
    });
    let data = [];
    for (let i = 0; i < squares.length; i++) {
      for (let j = 0; j < squares[i].length; j++) {
        if (squares[i][j] && squares[i][j].ids.length > 0) {
          data.push(squares[i][j]);
        }
      }
    }
    return { data };
  }
);

export const getOccupiedSquareGrid = createSelector(
  getSquareGrid,
  getAllScatterPlotDataBySquareBin,
  getZReducer,
  (grid, binned, reducer) => {
    if (!binned) return undefined;
    let data = [];
    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;
    binned.data.forEach((d) => {
      let index = grid.data.findIndex((o) => o.id === d.id);
      let cell = {
        id: d.id,
        i: d.i,
        j: d.j,
        x: grid.data[index].x,
        y: grid.data[index].y,
        width: grid.width,
        height: grid.height,
        ids: d.ids,
      };
      data.push(cell);
      let len = d.zs.length;
      for (let i = 0; i < len; i++) {
        if (d.zs[i] < min) min = d.zs[i];
      }
      let z = reducer.func(d.zs);
      max = Math.max(max, z);
      min = Math.min(min, z);
    });
    let range = [min, max];
    let newGrid = immutableUpdate(grid, { data, range });
    return newGrid;
  }
);

export const getScatterPlotDataBySquareBin = createSelector(
  getOccupiedSquareGrid,
  getScatterPlotData,
  getZReducer,
  getZScale,
  (grid, scatterData, reducer, scale) => {
    if (!grid) return undefined;
    let zScale = d3[scale]().domain(grid.range).range([0, grid.width]);
    let squares = [];
    grid.data.forEach((d) => {
      squares[d.i] = squares[d.i] || [];
      squares[d.i][d.j] = {
        id: d.id,
        i: d.i,
        j: d.j,
        ids: [],
        zs: [],
        indices: [],
      };
    });
    scatterData.data.forEach((datum) => {
      let coords = setCoords(datum, grid);
      if (squares[coords[0]] && squares[coords[0]][coords[1]]) {
        squares[coords[0]][coords[1]].indices.push(datum.index);
        squares[coords[0]][coords[1]].ids.push(datum.id);
        squares[coords[0]][coords[1]].zs.push(datum.z);
      }
    });
    let data = [];
    for (let i = 0; i < squares.length; i++) {
      if (squares[i]) {
        for (let j = 0; j < squares[i].length; j++) {
          if (squares[i][j] && squares[i][j].zs.length > 0) {
            let z = zScale(reducer.func(squares[i][j].zs));
            z = z > 1 ? z : 1;
            let offset = (grid.width - z) / 2;
            let index = grid.data.findIndex((o) => o.id === squares[i][j].id);
            squares[i][j].x = grid.data[index].x + offset;
            squares[i][j].y = grid.data[index].y + offset;
            squares[i][j].width = grid.width - offset * 2;
            squares[i][j].height = squares[i][j].width;
            squares[i][j].index = 0;
            data.push(squares[i][j]);
          }
        }
      }
    }
    return { data, grid, squares };
  }
);

export const getSelectedSquareGrid = createSelector(
  getOccupiedSquareGrid,
  getSelectedRecordsAsObject,
  getScatterPlotDataBySquareBin,
  (grid, records, scatterData) => {
    if (!scatterData) return undefined;
    let data = [];
    grid.data.forEach((cell) => {
      let newCell = Object.assign({}, cell);
      let index = scatterData.data.findIndex((o) => o.id === cell.id);
      if (index != -1) {
        newCell.ids = scatterData.data[index].ids;
        let count = 0;
        let len = newCell.ids.length;
        let i = 0;
        for (let i = 0; i < len; i++) {
          if (records[newCell.ids[i]]) count++;
        }
        newCell.selected = count;
      } else {
        newCell.ids = [];
      }
      data.push(newCell);
    });
    let newGrid = immutableUpdate(grid, { data });
    return newGrid;
  }
);

export const getScatterPlotDataBySquareBinByCategory = createSelector(
  getOccupiedSquareGrid,
  getScatterPlotDataBySquareBin,
  getBinsForCat,
  getFilteredDataForCat,
  getColorPalette,
  getZReducer,
  getZScale,
  getPlotScale,
  (
    grid = {},
    scatterData,
    bins,
    categories,
    palette,
    reducer,
    scale,
    plotScale
  ) => {
    if (!scatterData) return undefined;
    let zScale = d3[scale]().domain(grid.range).range([0, grid.width]);
    let keys = {};
    let data = [];
    if (bins.length == 0) return scatterData;
    bins.forEach((bin, i) => {
      bin.keys.forEach((key) => {
        keys[key] = i;
      });
    });
    let byCat = [];
    let max = Number.NEGATIVE_INFINITY;
    let min = Number.POSITIVE_INFINITY;
    let squareArrays = [];
    scatterData.data.forEach((square) => {
      let squareData = [];
      bins.forEach((bin, i) => {
        squareData[i] = {
          id: square.id + "_" + i,
          cellId: square.id,
          index: i,
          i: square.i,
          j: square.j,
          color: palette.colors[i],
          ids: [],
          zs: [],
          indices: [],
        };
      });
      square.indices.forEach((index, i) => {
        let currentCell = squareData[keys[categories.values[index]]];
        // if (currentCell){
        currentCell.ids.push(square.ids[i]);
        currentCell.indices.push(index);
        currentCell.zs.push(square.zs[i]);
        // }
      });
      let squareArray = squareData.filter((obj) => obj.ids.length > 0);
      squareArray.forEach((s) => {
        let value = reducer.func(s.zs);
        max = Math.max(max, value);
        min = Math.min(min, value);
      });
      squareArrays.push(squareArray);
    });

    zScale.domain([min, max]);
    squareArrays.forEach((squareArray) => {
      squareArray.forEach((s) => {
        s.z = zScale(reducer.func(s.zs)) * plotScale;
        s.z = s.z > plotScale ? s.z : plotScale;
        let offset = (grid.width - s.z) / 2;
        let index = grid.data.findIndex((o) => o.id === s.cellId);
        s.x = grid.data[index].x + offset;
        s.y = grid.data[index].y + offset;
        s.width = grid.width - offset * 2;
        s.height = s.width;
      });
      data = data.concat(squareArray.sort((a, b) => b.z - a.z));
    });
    return { data, grid, zScale };
  }
);

export const getSquareGridScale = createSelector(
  getScatterPlotDataBySquareBinByCategory,
  (binnedData) => binnedData.zScale
);

export const getBinnedDataByCategoryByAxis = createSelector(
  getScatterPlotDataBySquareBinByCategory,
  getBinsForCat,
  (binnedData, bins) => {
    if (!binnedData) return undefined;
    let iBinned = [{}];
    let jBinned = [{}];
    for (let index = 0; index < bins.length; index++) {
      iBinned[index] = {};
      jBinned[index] = {};
    }
    binnedData.data.forEach((d) => {
      if (!iBinned[d.index][d.i]) iBinned[d.index][d.i] = { zs: [], ids: [] };
      if (!jBinned[d.index][d.j]) jBinned[d.index][d.j] = { zs: [], ids: [] };
      iBinned[d.index][d.i].zs = iBinned[d.index][d.i].zs.concat(d.zs);
      //iBinned[d.index][d.i].ids = iBinned[d.index][d.i].ids.concat(d.ids)
      //iBinned[d.index][d.i].indices = iBinned[d.index][d.i].indices.concat(d.indices)
      jBinned[d.index][d.j].zs = jBinned[d.index][d.j].zs.concat(d.zs);
      //jBinned[d.index][d.j].ids = jBinned[d.index][d.j].ids.concat(d.ids)
      //jBinned[d.index][d.j].indices = jBinned[d.index][d.j].indices.concat(d.indices)
    });
    let data = [];
    for (let index = 0; index < Math.max(bins.length, 1); index++) {
      Object.keys(iBinned[index]).forEach((i) => {
        let obj = { axis: "x", bin: i, index: index };
        obj.zs = iBinned[index][i].zs;
        //obj.ids = iBinned[index][i].ids
        //obj.indices = iBinned[index][i].indices
        data.push(obj);
      });
      Object.keys(jBinned[index]).forEach((j) => {
        let obj = { axis: "y", bin: j, index: index };
        obj.zs = jBinned[index][j].zs;
        //obj.ids = jBinned[index][j].ids
        //obj.indices = jBinned[index][j].indices
        data.push(obj);
      });
    }
    let grid = binnedData.grid;
    return { data, grid };
  }
);

const getAxis = (state, axis) => axis;
export const getBinnedLinesByCategoryForAxis = createSelector(
  getAxis,
  getBinnedDataByCategoryByAxis,
  getFilteredDataForCat,
  getColorPalette,
  getSideMax,
  getZReducer,
  getZScale,
  (axis, binnedData, categories, palette, side, reducer, scale) => {
    if (!binnedData) return undefined;
    let min = 0;
    let max = side;
    if (side && side > 0) {
      max = side;
    } else {
      side = Number.NEGATIVE_INFINITY;
    }
    binnedData.data.forEach((d) => {
      let z = reducer.func(d.zs);
      max = Math.max(max, z);
    });

    let zScale = d3.scaleLinear().domain([0, max]).range([0, 300]);
    // zScale = d3.scaleSqrt().clamp(1).domain([1,max]).range([0,300])
    let paths = [{ name: "all", color: "#999999", path: "M0 300" }];
    if (categories.values.length > 0) {
      for (let i = 0; i < categories.keys.length; i++) {
        paths[i] = {
          name: categories.keys[i],
          color: palette.colors[i],
          path: "M0 300",
        };
      }
    }
    let zs = {};
    let xs = {};
    let bins = {};
    binnedData.data.forEach((d) => {
      if (d.axis == axis) {
        bins[d.index] = bins[d.index] || 0;
        zs[d.index] = zs[d.index] ? zs[d.index] : 300;
        let x0 = d.bin * binnedData.grid.width + 50;
        let x1 = x0 + binnedData.grid.width;
        let z = 300 - zScale(reducer.func(d.zs));
        if (bins[d.index] == 0) {
          paths[d.index].path += " L" + x0 + " " + 300;
        } else if (bins[d.index] != d.bin) {
          paths[d.index].path += " L" + xs[d.index] + " " + 300;
          paths[d.index].path += " L" + x0 + " " + 300;
        }
        paths[d.index].path += " L" + x0 + " " + z;
        paths[d.index].path += " L" + x1 + " " + z;
        zs[d.index] = z;
        xs[d.index] = x1;
        bins[d.index] = d.bin * 1 + 1;
      }
    });
    Object.keys(xs).forEach((index) => {
      //if (xs[index] < 1000){
      paths[index].path += " L" + xs[index] + " " + 300;
      paths[index].path += " L1000 " + 300;
      //}
    });
    return { paths, zScale };
  }
);
