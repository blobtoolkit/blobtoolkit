import { createAction, handleAction, handleActions } from 'redux-actions'
import { createSelector } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import { getIdentifiers, fetchIdentifiers } from './identifiers'
import immutableUpdate from 'immutable-update';
import deep from 'deep-get-set'
import shallow from 'shallowequal'
import store from '../store'
import qs from 'qs'
import { filterToList, getFilteredList, getUnfilteredList } from './filter'
import { fetchRawData, getAllActiveFields } from './field'
import queryToStore from '../querySync'
import { selectNone, addRecords, setSelectSource, changeSelectPolygon } from './select'
import { getQueryString, setHashString } from './location'
import { getSelectedDatasetMeta } from './dataset'

export const addList = createAction('ADD_LIST')
export const editList = createAction('EDIT_LIST')

export const lists = handleActions(
  {
    ADD_LIST: (state, action) => (
      immutableUpdate(state, {
        byId: { [action.payload.id]: action.payload },
        allIds: [...state.allIds, action.payload.id]
      })
    ),
    EDIT_LIST: (state, action) => {
      let id = action.payload.id
      let fields = Object.keys(action.payload).filter((key)=>{return key != 'id'})
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

const getListOfLists = state => state.lists

const arraysEqual = (a,b) => {
  if (a === b) return true;
  if (a == null || b == null) return false;
  if (a.length != b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

export const getLists = createSelector(
  getListOfLists,
  getFilteredList,
  getQueryString,
  getUnfilteredList,
  (lol,list,qStr,all) => {
    let ret = [
      {id:'current',list,params:list.params = qs.parse(qStr)}
    ]
    if (!arraysEqual(list,all)){
      ret.unshift({id:'all',list:all,params:{}})
    }
    ret = ret.concat(lol.allIds.map(id => lol.byId[id]))
    return ret
  }
)

//http://localhost:8080/view/Arthropoda/dataset/ACVV01/blob?
//gc--Min=0.368&
//SRR026696_cov--Min=0.01&
//bestsumorder_phylum--Keys=4%2C6%2C3%2C7%2C13%2C5%2C8%2C9%2C11%2C10%2C12%2C16%2C17%2C14%2C15%2C18%2C19%2C20%2C21%2C0&
//length--Min=1850&
//color1=rgba%280%2C0%2C0%2C1%29&
//color2=rgba%28248%2C231%2C28%2C1%29&
//selection--Active=true#Lists
export const updateSelectedList = createAction('UPDATE_SELECTED_LIST')
export const chooseList = (id,select) => {
  return async function (dispatch) {
    let state = store.getState()
    let list = getListById(state,id)
    if (list.hasOwnProperty('polygon')){
      await dispatch(changeSelectPolygon(list.polygon))
      if (list.hasOwnProperty('identifiers')){
        await dispatch(setSelectSource('circle'))
      }
      else {
        await dispatch(setSelectSource('list'))
      }

    }
    await dispatch(selectNone())
    let values = Object.assign({},list.params)
    await dispatch(queryToStore({values,searchReplace:true}))
    let fields = getAllActiveFields(store.getState())
    Object.keys(fields).forEach(async (field) =>{
      await dispatch(fetchRawData(field))
    })
    if (select){
      await dispatch(addRecords(list.list))
    }
    dispatch(setHashString('Filters'))


  }
}
export const selectedList = handleAction(
  'UPDATE_SELECTED_LIST',
  (state,action) => (
    action.payload
  ),
  null
)
const createSelectorForSelectedList = byIdSelectorCreator();
export const getSelectedList = state => state.selectedList


const createSelectorForListId = byIdSelectorCreator();
const _getListIdAsMemoKey = (state, id) => id;
const getList = (state, id) => state.lists ? (state.lists.byId[id] || {id}) : {id};

export const getListById = createSelectorForListId(
  _getListIdAsMemoKey,
  getList,
  getFilteredList,
  getUnfilteredList,
  getQueryString,
  (list,filtered,all,qStr) => {
    if (!list.list){
      if (list.id == 'all'){
        list.list = all
        list.params = {}
      }
      else {
        list.list = filtered
        list.params = qs.parse(qStr)
      }
    }
    return list
  }
);
// const getListById = (state,id) => {
//   let list = state.lists.byId[id]
//   console.log(id)
//   list = list || state.lists.byId[getSelectedList(state)]
//   if (!list){
//     if (id == 'all')
//       list = {id,list:getFilteredList(state)}
//     list.params = qs.parse(getQueryString(state))
//   }
//   return list
// }
const createSelectorForListIdentifiers = byIdSelectorCreator();
export const getIdentifiersForList = createSelectorForListIdentifiers(
  _getListIdAsMemoKey,
  getListById,
  getIdentifiers,
  (list,ids) => {
    let ret = [];
    (list.list || []).forEach(index => {ret.push(ids[index])})
    return ret;
  }
)

export function createList(val) {
  return function(dispatch){
    val.id = val.id.replace(/\s/g,'_')
    let list = val
    let state = store.getState()
    list.params = qs.parse(getQueryString(state))
    dispatch(addList(list))
  }
}

export function uploadedFileToList(acceptedFiles) {
  return function(dispatch){
    acceptedFiles.forEach(file => {
      const reader = new FileReader();
      reader.onload = () => {
          const fileAsBinaryString = reader.result;
          let obj = JSON.parse(fileAsBinaryString)
          dispatch(fetchIdentifiers()).then(ids=>{
            obj.list = []
            obj.id = 'user_'+obj.id
            if (ids.payload) ids = ids.payload
            let set = {}
            if (obj.hasOwnProperty('identifiers')){
              obj.identifiers.forEach(id=>{
                set[id]=true
              })
              for (let i = 0,id; id = ids[i]; i++){
                if (set[id]){
                  obj.list.push(i)
                }
              }
              delete obj.identifiers
            }
            else {
              obj.list = ids
              obj.params['selection--Active'] = false
            }
            dispatch(addList(obj))
          })
      };
      reader.onabort = () => console.log('file reading was aborted');
      reader.onerror = () => console.log('file reading has failed');
      reader.readAsBinaryString(file);
      window.URL.revokeObjectURL(file.preview);
    });
  }
}

export const listReducers = {
  lists,
  selectedList
}
