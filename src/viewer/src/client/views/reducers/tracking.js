import { createAction, handleAction, handleActions } from "redux-actions";

import ReactGA from "react-ga";
import { createSelector } from "reselect";
import store from "../store";

const getGaId = () => {
  if (
    window &&
    window.process &&
    window.process.ENV &&
    window.process.ENV.BTK_GA_ID
  ) {
    return window.process.ENV.BTK_GA_ID;
  }
  return GA_ID || "UA-000000-01";
};

const trackingId = GA_ID || "UA-000000-01";

const GA_CONFIG = {
  trackingId,
  debug: true,
  gaOptions: {
    cookieDomain: "none",
  },
};

export const setCookieConsent = createAction("SET_COOKIE_CONSENT");

export const cookieConsent = handleAction(
  "SET_COOKIE_CONSENT",
  (state, action) => action.payload,
  false
);

export const getCookieConsent = (state) => state.cookieConsent;

export const setAnalytics = createAction("SET_ANALYTICS");

export const analytics = handleAction(
  "SET_ANALYTICS",
  (state, action) => action.payload,
  false
);

export const getAnalytics = (state) => state.analytics;

export const trackPage = (page) => {
  ReactGA.set({
    page,
  });
  ReactGA.pageview(page);
};

export const startAnalytics = () => {
  return async function (dispatch) {
    let state = store.getState();
    let consent = getCookieConsent(state);
    let analytics = getAnalytics(state);
    if (consent && !analytics) {
      if (trackingId && trackingId != "UA-000000-01") {
        dispatch(setAnalytics(true));
        ReactGA.initialize(GA_CONFIG);
        trackPage(window.location.pathname + window.location.search);
      }
    }
  };
};

export const trackingReducers = {
  cookieConsent,
  analytics,
};
