import { createAction, handleAction, handleActions } from 'redux-actions'
import { createSelector, createSelectorCreator } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import immutableUpdate from 'immutable-update';
import deep from 'deep-get-set'
import shallow from 'shallowequal'
import store from '../store'
import { queryToStore, colorToRGB, qsDefault, userColors } from '../querySync'
import { getQueryValue } from '../reducers/location'

export const addPalette = createAction('ADD_PALETTE')
export const editPalette = createAction('EDIT_PALETTE')

const qsPalette = () => {
  let colors = userColors.slice(0)
  for (let i = 0; i < colors.length; i++){
    let qsColor = qsDefault('color'+i)
    if (qsColor){
      colors[i] = colorToRGB(qsColor) || colors[i]
    }
  }
  return colors
}

export const palettes = handleActions(
  {
    ADD_PALETTE: (state, action) => (
      immutableUpdate(state, {
        byId: { [action.payload.id]: action.payload },
        allIds: [...state.allIds, action.payload.id]
      })
    ),
    EDIT_PALETTE: (state, action) => {
      let id = action.payload.id
      let arr = action.payload[id].slice(0)
      state.byId[id].forEach((col,i)=>{
        if (!arr[i]) arr[i] = col
      })
      return immutableUpdate(state, {
        byId: {
          [id]: arr
        }
      })
    }
  },
  {
    byId: {
      default:['rgb(166,206,227)','rgb(31,120,180)','rgb(178,223,138)','rgb(51,160,44)','rgb(251,154,153)','rgb(227,26,28)','rgb(253,191,111)','rgb(255,127,0)','rgb(202,178,214)','rgb(106,61,154)'],
      user: qsPalette()
    },
    allIds: ['default','user']
  }
)

export const getSelectedPalette = state => state.selectedPalette

export const selectPalette = createAction('SELECT_PALETTE')
export const selectedPalette = handleAction(
  'SELECT_PALETTE',
  (state, action) => (
    action.payload
  ),
  qsDefault('palette')
)

export const choosePalette = (palette) => {
  return function (dispatch) {
    let values = {palette}
    dispatch(queryToStore({values}))
    //dispatch(selectPalette(palette))
  }
}

export const chooseColors = (colors) => {
  return function (dispatch) {
    let existing = getColorPalette(store.getState())
    let values = {colors:{colors,existing}}
    dispatch(queryToStore({values}))
    //dispatch(selectPalette(palette))
  }
}

export const getAllPalettes = state => state.palettes

export const getColorPalette = createSelector(
  getSelectedPalette,
  getAllPalettes,
  (id, palettes, qsPalette) => {
    id = qsPalette || id
    let colors = palettes ? palettes.byId[id] : []
    return {id,colors}
  }
)

export const getUserPalette = createSelector(
  getAllPalettes,
  (palettes) => {
    let id = 'user'
    let colors = palettes ? palettes.byId[id] : []
    return {id,colors}
  }
)

export const getDefaultPalette = createSelector(
  getAllPalettes,
  (palettes) => {
    let id = 'default'
    let colors = palettes ? palettes.byId[id] : []
    return {id,colors}
  }
)


/* Coolors Exported Palette - coolors.co/d7cdcc-ffffff-59656f-9c528b-1d1e2c */
/* RGB */
export const schemeColors = {
  darkColor: 'rgb(49, 50, 63)',
  lightColor: 'rgb(255, 255, 255)',
  shadeColor: 'rgb(89, 101, 111)',
  deepColor: 'rgb(65,74,81)',
  highlightColor: 'rgb(156, 82, 139)',
  halfHighlightColor: 'rgba(156, 82, 139, 0.5)',
  paleColor: 'rgb(215, 205, 204)',
  brightColor: 'rgb(255,255,30)',
  clearColor: 'rgba(255,255,255,0)'
}
export const setColorScheme = createAction('SET_COLOR_SCHEME')
export const colorScheme = handleAction(
  'SET_COLOR_SCHEME',
  (state, action) => (
    action.payload
  ),
  schemeColors
)

export const getColorScheme = state => state.colorScheme


export const colorReducers = {
  palettes,
  selectedPalette,
  colorScheme
}
