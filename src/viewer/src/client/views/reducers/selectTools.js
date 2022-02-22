import { createAction, handleAction, handleActions } from 'redux-actions'
import { createSelector } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import { getFilteredList } from './filter'
import { getCurrentDatasetMeta } from './dataset'
import immutableUpdate from 'immutable-update';
import deep from 'deep-get-set'
import store from '../store'
import { addRecords, selectNone, getSelectedRecordsAsObject } from './select'


export function selectAll() {
  return dispatch => {
    let state = store.getState()
    let meta = getCurrentDatasetMeta(state)
    let all = getFilteredList(state)
    dispatch(addRecords(all))
  }
}

export function invertSelection() {
  return dispatch => {
    let state = store.getState()
    let meta = getCurrentDatasetMeta(state)
    // let all = getFilteredList(state)
    let current = getSelectedRecordsAsObject(state)
    let all = []
    let rev = []
    for (let i = 0; i < meta.records; i++){
      if (!current[i]){
        rev.push(i)
      }
      else {
        all.push(i)
      }
    }
    dispatch(selectNone())
    dispatch(addRecords(rev))
  }
}

export function selectInverse() {
  return dispatch => {
    let state = store.getState()
    let meta = getCurrentDatasetMeta(state)
    let all = []
    for (let i = 0; i < meta.records; i++){
      all.push(i)
    }
    dispatch(addRecords(all))
  }
}
