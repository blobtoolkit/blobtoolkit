export const getApiURL = () => {
  if (
    window &&
    window.process &&
    window.process.ENV &&
    window.process.ENV.BTK_API_URL
  ) {
    return window.process.ENV.BTK_API_URL;
  }
  return API_URL || "api/v1";
};

export const apiUrl = getApiURL();
