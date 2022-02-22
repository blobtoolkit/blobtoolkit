import { createAction, handleAction, handleActions } from 'redux-actions'
import { createSelector } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import immutableUpdate from 'immutable-update';
import deep from 'deep-get-set'
import store from '../store'

export const setDimension = createAction('SET_DIMENSION')

export const dimensions = handleAction(
  'SET_DIMENSION',
  (state, action) => {
    let id = action.payload.id
    let fields = Object.keys(action.payload).filter((key)=>{return key != 'id'})
    return immutableUpdate(state, {
      byId: {
        [id]: Object.assign(...fields.map(f => ({[f]: action.payload[f]})))
      }
    })
  },
  {
    byId: {
      preview:{
        width:400,
        height:400/3.75,
        xDomain:[0,400],
        yDomain:[0,400/3.75],
        xScale:'scaleLinear',
        yScale:'scaleLinear'
      }
    }
  }
)

const createSelectorForDimensionId = byIdSelectorCreator();

const _getDimensionIdAsMemoKey = (state, dimensionId) => dimensionId;
const getDimensions = (state, dimensionId) => state.dimensions.byId[dimensionId];

export const getDimensionsbyDimensionId = createSelectorForDimensionId(
  _getDimensionIdAsMemoKey,
  getDimensions,
  (dimensions) => dimensions
);

export const getPreviewDimensions = createSelector(
  (state) => getDimensionsbyDimensionId(state,'preview'),
  (dimensions) => dimensions
)

export const dimensionReducers = {
  dimensions
}
