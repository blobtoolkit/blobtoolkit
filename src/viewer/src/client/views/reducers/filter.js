import { createAction, handleAction, handleActions } from 'redux-actions'
import { createSelector } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import immutableUpdate from 'immutable-update';
import deep from 'deep-get-set'
import shallow from 'shallowequal'
import store from '../store'
import { getQueryValue } from './location'
import { queryToStore } from '../querySync'
import qs from 'qs'
import { getSelectedRecords } from './select'
import { getDatasetIsActive } from './repository'
import { getSelectedDatasetMeta } from './dataset'

export const addFilter = createAction('ADD_FILTER')
export const editFilter = createAction('EDIT_FILTER')

function isNumeric(n) {
  if ((typeof n === 'undefined') || n == 'NaN') return false
  return !isNaN(parseFloat(n)) && isFinite(n)
}

export const filters = handleActions(
  {
    ADD_FILTER: (state, action) => (
      immutableUpdate(state, {
        byId: { [action.payload.id]: action.payload },
        allIds: [...state.allIds, action.payload.id]
      })
    ),
    EDIT_FILTER: (state, action) => {
      let id = action.payload.id
      let fields = Object.keys(action.payload).filter((key)=>{return key != 'id'})
      let limit = []
      if (state.byId[id] && action.payload.limit){
        if (action.payload.limit.length == 1){
          limit[0] = isNumeric(action.payload.limit[0]) ? action.payload.limit[0] : Number.NEGATIVE_INFINITY
          limit[1] = state.byId[id] ? state.byId[id].limit[1] : Number.POSITIVE_INFINITY
        }
        else if (typeof action.payload.limit[0] == 'undefined'){
          limit[0] = state.byId[id] ? state.byId[id].limit[0] : Number.NEGATIVE_INFINITY
          limit[1] = isNumeric(action.payload.limit[1]) ? action.payload.limit[1] : Number.POSITIVE_INFINITY
        }
        else {
          limit = action.payload.limit.slice()
        }
        action.payload.limit = limit
      }
      if (state.byId[id] && action.payload.range){
        let range = []
        if (action.payload.range.length == 1){
          range[0] = isNumeric(action.payload.range[0]) ? action.payload.range[0] : Number.NEGATIVE_INFINITY
          range[1] = state.byId[id] ? state.byId[id].range[1] : Number.POSITIVE_INFINITY
        }
        else if (typeof action.payload.range[0] == 'undefined'){
          range[0] = state.byId[id] ? state.byId[id].range[0] : Number.NEGATIVE_INFINITY
          range[1] = isNumeric(action.payload.range[1]) ? action.payload.range[1] : Number.POSITIVE_INFINITY
        }
        else {
          range = action.payload.range.slice()
          if (isNaN(range[0])) range[0] = state.byId[id] ? state.byId[id].range[0] : Number.NEGATIVE_INFINITY
          if (isNaN(range[1])) range[1] = state.byId[id] ? state.byId[id].range[1] : Number.NEGATIVE_INFINITY
        }
        if (range.length > 0){
          if (limit.length == 0){
            limit = (state.byId[id] ? (state.byId[id].limit || state.byId[id].range) : range).slice()
            if (state.byId[id] && action.payload.limit){
              action.payload.limit = limit
            }
          }
          action.payload.range[0] = Math.max(range[0],limit[0])
          action.payload.range[1] = Math.min(range[1],limit[1])
        }
      }

      return immutableUpdate(state, {
        byId: {
          [id]: Object.assign(...fields.map(f => ({[f]: action.payload[f]})))
        }
      })
    }
  },
  {
    byId: {},
    allIds: []
  }
)

const filterRangeToList = (low,high,arr,list,invert) => {
  let ret = []
  let len = list.length
  for (var i = 0; i < len; i++){
    if (invert){
      if (arr[list[i]] < low || arr[list[i]] > high){
        ret.push(list[i]);
      }
    }
    else {
      if (arr[list[i]] >= low && arr[list[i]] <= high){
        ret.push(list[i]);
      }
    }
  }
  return ret
}

const filterCategoriesToList = (keys,arr,list,invert) => {
  let ret = []
  let len = list.length
  for (var i = 0; i < len; i++){
    if (invert){
      if (keys.includes(arr[list[i]])){
        ret.push(list[i]);
      }
    }
    else {
      if (!keys.includes(arr[list[i]])){
        ret.push(list[i]);
      }
    }
  }
  return ret
}

const filterArrayToList = (arr,list) => {
  let ret = []
  let a=0, l=0;
  while (a < arr.length && l < list.length){
    if (arr[a] < list[l] ){
      a++
    }
    else if (arr[a] > list[l]){
      l++
    }
    else {
      ret.push(arr[a])
      a++
      l++
    }
  }
  return ret
}

const filterArrayFromList = (arr,list) => {
  let ret = []
  let a=0, l=0;
  while (a < arr.length && l < list.length){
    if (arr[a] < list[l] ){
      a++
    }
    else if (arr[a] > list[l]){
      ret.push(list[l])
      l++
    }
    else {
      a++
      l++
    }
  }
  return ret
}

const keyed = (o,k) => Object.prototype.hasOwnProperty.call(o,k)

export function filterToList(readQueryString) {
  return function(dispatch){
    let state = store.getState();
    let filters = state.filters.byId;
    let fields = state.fields.byId;
    let data = state.rawData.byId;
    let parsed = qs.parse(state.queryString)
    let count = getSelectedDatasetMeta(state).records
    let list = fields['selection'].active ? state.selectedRecords : undefined
    let all = []
    if (!list || list.length == 0 || filters['selection'].invert){
      for (let i = 0; i < count; i++){
        all.push(i)
      }
    }
    if (!list || list.length == 0){
      list = all
    }
    else if (fields['selection'].active && filters['selection'].invert){
      list = filterArrayFromList(list,all)
    }
    let remove = []
    let values = {}
    state.filters.allIds.forEach(id => {
      if (data[id] && fields[id] && fields[id].active && filters[id]){
        if (filters[id].type == 'range'){
          let minstr = id+'--Min'
          let maxstr = id+'--Max'
          let invstr = id+'--Inv'
          let range = filters[id].range
          let limit = fields[id].range
          let qmin = 1 * getQueryValue(minstr)
          let qmax = 1 * getQueryValue(maxstr)
          if (!shallow(range,limit)){
            list = filterRangeToList(range[0],range[1],data[id].values,list,filters[id].invert)
            if (range[0] > limit[0]){
              values[minstr] = range[0]
            }
            else if (keyed(parsed,minstr)){
              remove.push(minstr)
            }
            if (range[1] < limit[1]){
              values[maxstr] = range[1]
            }
            else if (keyed(parsed,maxstr)){
              remove.push(maxstr)
            }
            if (filters[id].invert){
              values[invstr] = true
            }
            else if (getQueryValue(invstr)){
              values[invstr] = false
            }
          }
          else {
            if (keyed(parsed,minstr)){
              remove.push(minstr)
            }
            if (keyed(parsed,maxstr)){
              remove.push(maxstr)
            }
          }
        }
        else if (filters[id].type == 'list' || filters[id].type == 'category'){
          if (data[id]){
            let keys = []
            filters[id].keys.forEach(key=>{
              if (isNaN(key)){
                key = data[id].keys.indexOf(key)
              }
              if (key > -1){
                keys.push(key)
              }
            })
            list = filterCategoriesToList(keys,data[id].values,list,filters[id].invert)
            let keystr = id+'--Keys'
            let invstr = id+'--Inv'
            if (keys.length > 0){
              values[keystr] = keys.join()
              if (filters[id].invert){
                values[invstr] = true
              }
              else {
                values[invstr] = false
              }
            }
            else {
              if (keyed(parsed,keystr)){
                remove.push(keystr)
              }
              if (keyed(parsed,invstr)){
                remove.push(invstr)
              }
            }
          }
        }
      }
    })
    dispatch(queryToStore({values,remove}))
  }
}

const getAllFilters = createSelector(
  (state) => state.filters,
  filters => {
    return filters
  }
)

const createSelectorForFilterId = byIdSelectorCreator();

const _getFilterIdAsMemoKey = (state, filterId) => filterId;
export const getMetaDataForFilter = (state, filterId) => state.filters ? state.filters.byId[filterId] : {};


const getAllFields = createSelector(
  (state) => state.fields ? state.fields.byId : {},
  fields => {
    return fields
  }
)

const getAllRawData = createSelector(
  (state) => state.rawData ? state.rawData.byId : {},
  rawData => {
    return rawData
  }
)

export const getAllActiveFilters = createSelector(
  getAllFilters,
  getAllFields,
  getAllRawData,
  (filters, fields, data) => {
    let active = []
    filters.allIds.forEach(id=>{
      let filter = filters.byId[id]
      let field = fields[id]
      let raw = data[id]
      if (field.active){
        if (filter.type == 'range' && (filter.range[0] > field.range[0] || filter.range[1] < field.range[1])){
          if (filter.range[0] > field.range[0] && filter.range[1] < field.range[1]){
            if (filter.invert){
              active.push({id,type:'out',range:filter.range})
            }
            else {
              active.push({id,type:'in',range:filter.range})
            }
          }
          else if (filter.range[0] > field.range[0]){
            if (filter.invert){
              active.push({id,type:'lt',value:filter.range[0]})
            }
            else {
              active.push({id,type:'gt',value:filter.range[0]})
            }
          }
          else if (filter.range[1] < field.range[1]){
            if (filter.invert){
              active.push({id,type:'gt',value:filter.range[1]})
            }
            else {
              active.push({id,type:'lt',value:filter.range[1]})
            }
          }
        }
        else if (filter.type == 'category' && filter.keys.length > 0){
          let type = 'cat'
          if (filter.invert){
            type = 'nocat'
          }
          active.push({id,type,list:filter.keys.map(k=>raw.keys[k])})
        }
      }

    })
    return active
  }
)



export const getActiveSelection = createSelector(
  getAllFields,
  getSelectedRecords,
  (fields,list) => {
    if (fields['selection'] && fields['selection'].active){
      return list
    }
    return []
  }
)

export const getUnfilteredList = createSelector(
  getSelectedDatasetMeta,
  (meta) => {
    let all = []
    let count = meta.records || 0
    for (let i = 0; i < count; i++){
      all.push(i)
    }
    return all
  }
)


const getDatasetActive = createSelector(
  (state) => getDatasetIsActive(state),
  active => {
    return active
  }
)

export const getFilteredList = createSelector(
  getAllFilters,
  getAllFields,
  getAllRawData,
  getSelectedDatasetMeta,
  getActiveSelection,
  getDatasetActive,
  getUnfilteredList,
  (filters,fields,data,meta,list,active,all) => {
    let obj = {filters,fields,data,meta}
    let count = meta.records | 0
    if (!list || list.length == 0){
      list = all
    }
    else if (fields['selection'].active && filters.byId['selection'].invert){
      list = filterArrayFromList(list,all)
    }
    let values = {}
    if (filters){
      filters.allIds.forEach(id => {
        if (fields[id] && fields[id].active && filters.byId[id] && data[id]){
          if (filters.byId[id].type == 'range'){
            let range = filters.byId[id].range
            let limit = fields[id].range
            if (isNaN(range[0])) range[0] = limit[0]
            if (isNaN(range[1])) range[1] = limit[1]
            if (!shallow(range,limit)){
              list = filterRangeToList(range[0],range[1],data[id].values,list,filters.byId[id].invert)
            }
          }
          else if (filters.byId[id].type == 'list' || filters.byId[id].type == 'category'){
            if (data[id]){
              list = filterCategoriesToList(filters.byId[id].keys,data[id].values,list,filters.byId[id].invert)
            }
          }
        }
      })
    }
    return list
  }
)

export const filterReducers = {
  filters
}
