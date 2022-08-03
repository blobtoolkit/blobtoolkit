export const getGdprURL = () => {
  if (
    window &&
    window.process &&
    window.process.ENV &&
    window.process.ENV.BTK_GDPR_URL
  ) {
    return window.process.ENV.BTK_GDPR_URL;
  }
  return GDPR_URL;
};

export const gdprUrl = getGdprURL();
