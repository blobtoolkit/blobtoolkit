import { createSelector } from 'reselect'
import { byIdSelectorCreator } from './selectorCreators'
import { createBrowserHistory } from 'history'
import qs from 'qs'
const basename = BASENAME || ''
export const history = createBrowserHistory({
  basename
})

export default history;
