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

export const getTableFlag = () => {
  if (
    window &&
    window.process &&
    window.process.ENV &&
    window.process.ENV.BTK_DATASET_TABLE
  ) {
    return window.process.ENV.BTK_DATASET_TABLE ? true : false;
  }
  return DATASET_TABLE;
};

export const datasetTable = getTableFlag();

export const getMessage = () => {
  if (
    window &&
    window.process &&
    window.process.ENV &&
    window.process.ENV.BTK_MESSAGE
  ) {
    return window.process.ENV.BTK_MESSAGE
      ? JSON.stringify(window.process.ENV.BTK_MESSAGE)
      : false;
  }
  return MESSAGE;
};

export const message = getMessage();
