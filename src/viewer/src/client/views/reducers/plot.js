import { createAction, handleAction, handleActions } from "redux-actions";

import { byIdSelectorCreator } from "./selectorCreators";
import { createSelector } from "reselect";
import deep from "deep-get-set";
import { getDatasetMeta } from "../reducers/repository";
import { getQueryValue } from "./location";
import immutableUpdate from "immutable-update";
import { qsDefault } from "../querySync";
import store from "../store";

export const editPlot = createAction("EDIT_PLOT");
const defaultPlot = () => {
  return {};
};
export const plot = handleAction(
  "EDIT_PLOT",
  (state, action) => {
    if (!action.payload) return {};
    let fields = Object.keys(action.payload).filter((key) => {
      return key != "id";
    });
    if (!fields || fields.length == 0) return {};
    return immutableUpdate(
      state,
      Object.assign(...fields.map((f) => ({ [f]: action.payload[f] })))
    );
  },
  {}
);
export const getPlot = (state) => state.plot;

export const getMainPlot = createSelector(getPlot, (axes) => ({
  id: "default",
  axes,
}));

export const getXAxis = createSelector(
  (state) => getMainPlot(state).axes.x,
  (axis) => axis
);

export const getYAxis = createSelector(
  (state) => getMainPlot(state).axes.y,
  (axis) => axis
);

export const getZAxis = createSelector(
  (state) => getMainPlot(state).axes.z,
  (axis) => axis
);

export const getCatAxis = createSelector(
  (state) => getMainPlot(state).axes.cat,
  (axis) => axis
);

export const plotReducers = {
  plot,
};
