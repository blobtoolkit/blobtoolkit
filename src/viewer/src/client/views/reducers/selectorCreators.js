import { createSelectorCreator } from 'reselect'
import { handleAction } from 'redux-actions'
import shallow from 'shallowequal'
import immutableUpdate from 'immutable-update';
import store from '../store'

export const byIdSelectorCreator = () => {
  return createSelectorCreator((resultFunc) => {
    const memoAll = {};
    return (id, ...args) => {
      if (!memoAll[id]) {
        memoAll[id] = {};
      }
      const memo = memoAll[id];
      if (!shallow(memo.lastArgs, args)) {
        memo.lastArgs = args;
        memo.lastResult = resultFunc(...args);
      }
      return memo.lastResult;
    };
  })
}
