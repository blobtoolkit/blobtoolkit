import { createSelector } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import { getRawDataForFieldId, getDetailsForFieldId } from './field';
import { getBinsForFieldId, getBinsForCat } from './field';
import { getFilteredDataForFieldId } from './preview'
import { getMainPlot } from './plot';
import { getSelectedRecordsAsObject } from './select';
import { getMainPlotData, getScatterPlotData, getAllScatterPlotData, getFilteredDataForCat } from './plotData';
import { getZReducer, getZScale, getPlotResolution, getTransformFunction, getPlotScale } from './plotParameters';
import { getColorPalette } from './color';
import { drawPoints, pixel_to_oddr } from './hexFunctions'
import * as d3 from 'd3'
import immutableUpdate from 'immutable-update';

export const getHexGrid = createSelector(
  getPlotResolution,
  (res) => {
    let size = 900.1 // FIXME: magic number
    let width = size / res
    let jres = res * 3 / Math.sqrt(3)
    let radius = width / Math.sqrt(3)
    let height = radius * 2
    let data = [];
    for (let i = 0; i <= res; i++){
      for (let j = 0; j <= jres; j++){
        let points = drawPoints(i,j,radius,res)
        if (points) data.push({id:i+'_'+j,x:i,y:j,points})
      }
    }
    return { data,size,res,radius,height,width }
  }
)

export const getAllScatterPlotDataByHexBin = createSelector(
  getHexGrid,
  getAllScatterPlotData,
  (grid,scatterData) => {
    if (!scatterData) return undefined
    let hexes = []
    grid.data.forEach(d=>{
      hexes[d.x] = hexes[d.x] || []
      hexes[d.x][d.y] = {id:d.id,x:d.x,y:d.y,ids:[],zs:[],indices:[]}
    })
    scatterData.data.forEach(datum=>{
      let coords = pixel_to_oddr(datum.x,datum.y,grid.radius)
      if (hexes[coords.i] && hexes[coords.i][coords.j]){
        hexes[coords.i][coords.j].ids.push(datum.id)
        hexes[coords.i][coords.j].indices.push(datum.index)
        hexes[coords.i][coords.j].zs.push(datum.z)
      }
    })
    let data = [];
    for (let x = 0; x < hexes.length; x++){
      for (let y = 0; y < hexes[x].length; y++){
        if (hexes[x][y] && hexes[x][y].ids.length > 0){
          data.push(hexes[x][y])
        }
      }
    }
    return {data};
  }
)

export const getOccupiedHexGrid = createSelector(
  getHexGrid,
  getAllScatterPlotDataByHexBin,
  getZReducer,
  (grid,binned,reducer) => {
    if (!binned) return undefined
    let data = []
    let min = Number.POSITIVE_INFINITY
    let max = Number.NEGATIVE_INFINITY
    binned.data.forEach(d => {
      data.push({id:d.id,x:d.x,y:d.y,points:grid.data[grid.data.findIndex(o => o.id === d.id)].points})
      let len = d.zs.length
      for (let i = 0; i < len; i++){
        if (d.zs[i] < min) min = d.zs[i]
      }
      let z = reducer.func(d.zs)
      max = Math.max(max,z)
    })
    let range = [min,max]
    let newGrid = immutableUpdate(grid,{data,range})
    return newGrid;
  }
)


export const getScatterPlotDataByHexBin = createSelector(
  getOccupiedHexGrid,
  getScatterPlotData,
  getZReducer,
  getZScale,
  (grid,scatterData,reducer,scale) => {
    if (!grid || !scatterData) return undefined
    let zScale = d3[scale]().domain(grid.range).range([0,grid.radius])
    let hexes = []
    grid.data.forEach(d=>{
      hexes[d.x] = hexes[d.x] || []
      hexes[d.x][d.y] = {id:d.id,x:d.x,y:d.y,ids:[],zs:[],indices:[]}
    })
    scatterData.data.forEach(datum=>{
      let coords = pixel_to_oddr(datum.x,datum.y,grid.radius)
      if (hexes[coords.i] && hexes[coords.i][coords.j]){
        hexes[coords.i][coords.j].indices.push(datum.index)
        hexes[coords.i][coords.j].ids.push(datum.id)
        hexes[coords.i][coords.j].zs.push(datum.z)
      }
    })
    let data = [];
    for (let x = 0; x < hexes.length; x++){
      if (hexes[x]){
        for (let y = 0; y < hexes[x].length; y++){
          if (hexes[x][y] && hexes[x][y].ids.length > 0){
            let z = zScale(reducer.func(hexes[x][y].zs))
            z = z > 1 ? z : 1;
            hexes[x][y].points = drawPoints(x,y,grid.radius,grid.res,z/grid.radius)
            hexes[x][y].index = 0
            data.push(hexes[x][y])
          }
        }
      }
    }
    return {data,hexes,radius:grid.radius,grid};
  }
)

export const getSelectedHexGrid = createSelector(
  getOccupiedHexGrid,
  getSelectedRecordsAsObject,
  getScatterPlotDataByHexBin,
  (grid,records,scatterData) => {
    let data = []
    grid.data.forEach(cell => {
      let newCell = Object.assign({},cell)
      let index = scatterData.data.findIndex(o => o.id === cell.id)
      if (index != -1){
        newCell.ids = scatterData.data[index].ids
        let count = 0;
        let len = newCell.ids.length
        let i = 0
        for (let i = 0; i < len; i++){
          if (records[newCell.ids[i]]) count++
        }
        newCell.selected = count
      }
      else {
        newCell.ids = []
      }
      data.push(newCell)
    })
    let newGrid = immutableUpdate(grid,{data})
    return newGrid;
  }
)

export const getScatterPlotDataByHexBinByCategory = createSelector(
  getOccupiedHexGrid,
  getScatterPlotDataByHexBin,
  getBinsForCat,
  getFilteredDataForCat,
  getColorPalette,
  getZReducer,
  getZScale,
  getPlotScale,
  (grid,scatterData,bins,categories,palette,reducer,scale,plotScale) => {
    let zScale = d3[scale]().domain(grid.range).range([0,grid.radius])
    let keys = {}
    let data = []
    if (bins.length == 0) return scatterData
    bins.forEach((bin,i)=>{
      bin.keys.forEach(key=>{
        keys[key] = i
      })
    })
    let byCat = []
    let max = Number.NEGATIVE_INFINITY
    let min = Number.POSITIVE_INFINITY

    let hexArrays = []
    scatterData.data.forEach((hex)=>{
      let hexData = []

      bins.forEach((bin,i)=>{
        hexData[i] = {
          id:hex.id+'_'+i,
          x:hex.x,
          y:hex.y,
          color:palette.colors[i],
          ids:[],
          zs:[],
          indices:[]
        }
      })
      hex.indices.forEach((index,i) => {
        let currentCell = hexData[keys[categories.values[index]]]
        currentCell.ids.push(hex.ids[i])
        currentCell.indices.push(index)
        currentCell.zs.push(hex.zs[i])
      })
      let hexArray = hexData.filter(obj => obj.ids.length > 0);
      hexArray.forEach(h=>{
        let value = reducer.func(h.zs)
        max = Math.max(max,value)
        min = Math.min(min,value)
      })
      hexArrays.push(hexArray)

    })

    zScale.domain([min,max])
    hexArrays.forEach((hexArray)=>{
      hexArray.forEach(h=>{
        h.z = zScale(reducer.func(h.zs)) * plotScale
        h.z = h.z > plotScale ? h.z : plotScale;
        h.points = drawPoints(h.x,h.y,grid.radius,grid.res,h.z/grid.radius)
      })
      data = data.concat(hexArray.sort((a,b)=>b.z - a.z))

    })
    return {data, grid, zScale};
  }
)

export const getHexGridScale = createSelector(
  getScatterPlotDataByHexBinByCategory,
  (binnedData) => binnedData.zScale
)
