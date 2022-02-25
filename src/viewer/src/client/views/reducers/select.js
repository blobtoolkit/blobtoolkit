import { createAction, handleAction, handleActions } from 'redux-actions'
import { createSelector } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import { getIdentifiers } from './identifiers';
import immutableUpdate from 'immutable-update';
import deep from 'deep-get-set'
import store from '../store'


export const toggleSelection = createAction('TOGGLE_SELECTION')
const selectionDisplay = handleAction(
  'TOGGLE_SELECTION',
  (state, action) => (
    !action.payload
  ),
  true
)

const createSelectorForSelectionDisplay = byIdSelectorCreator();
export const getSelectionDisplay = state => state.selectionDisplay


export const addRecords = createAction('ADD_RECORDS')
export const removeRecords = createAction('REMOVE_RECORDS')
export const replaceRecords = createAction('REPLACE_RECORDS')
export const selectNone = createAction('SELECT_NONE')

export const selectedRecords = handleActions(
  {
    ADD_RECORDS: (state, action) => {
      let current = state || []
      let update = action.payload.slice(0)
      update.sort((a,b)=>a-b)
      let iLen = current.length
      let jLen = update.length
      let combined = []
      let i = 0
      let j = 0
      while (i < iLen && j < jLen) {
        if (current[i] < update[j]){
          combined.push(current[i])
          i++
        }
        else if (update[j] < current[i]){
          combined.push(update[j])
          j++
        }
        else {
          combined.push(current[i])
          i++
          j++
        }
      }
      while (i < iLen){
        combined.push(current[i])
        i++
      }
      while (j < jLen){
        combined.push(update[j])
        j++
      }
      return  combined
    },
    REMOVE_RECORDS: (state, action) => {
      let current = state
      let update = action.payload.slice(0)
      update.sort((a,b)=>a-b)
      let arr = []
      let iLen = current.length
      let jLen = update.length
      let i = 0
      let j = 0
      while (i < iLen && j < jLen) {
        if (current[i] < update[j]){
          arr.push(current[i])
          i++
        }
        else if (update[j] < current[i]){
          j++
        }
        else {
          i++
          j++
        }
      }
      while (i < iLen){
        arr.push(current[i])
        i++
      }
      return arr
    },
    REPLACE_RECORDS: (state, action) => {
      return action.payload.sort((a,b)=>a-b)
    },
    SELECT_NONE: (state, action) => {
      return []
    }
  },
  []
)

export function recordsByName(names) {
  return function (dispatch) {
    let state = store.getState();
    let existing = getIdentifiers(state)
    let records = []
    names.forEach(name=>{
      let index = existing.indexOf(name)
      if (index >= 0){
        records.push(index)
      }
    })
    dispatch(replaceRecords(records, true))
  }
}


export const getSelectedRecords = state => state.selectedRecords

export const getSelectedRecordsAsObject = createSelector(
  getSelectedRecords,
  getSelectionDisplay,
  (arr,display) => {
    let obj = {}
    if (!display) return obj
    arr.forEach(id=>{
      obj[id] = 1
    })
    return obj
  }
)

export const changeSelectPolygon = createAction('CHANGE_SELECT_POLYGON')
export const selectPolygon = handleAction(
  'CHANGE_SELECT_POLYGON',
  (state, action) => (
    action.payload
  ),
  []
)
export const getSelectPolygon = state => state.selectPolygon

export const setSelectSource = createAction('SET_SELECT_SOURCE')
export const selectSource = handleAction(
  'SET_SELECT_SOURCE',
  (state, action) => (
    action.payload
  ),
  'unknown'
)
export const getSelectSource = state => state.selectSource



export const selectReducers = {
  selectedRecords,
  selectionDisplay,
  selectPolygon,
  selectSource
}
