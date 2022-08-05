require("winston-daily-rotate-file");

const winston = require("winston");
const { createLogger, format, transports } = winston;

const config = require("../config/main");

const { combine, timestamp, prettyPrint, colorize, errors, printf } = format;

const errorTransport = new transports.DailyRotateFile({
  filename: config.errorLog,
  level: "error",
  format: combine(
    errors({ stack: true }),
    colorize(),
    timestamp(),
    prettyPrint()
  ),
  datePattern: "YYYY-MM-DD",
  zippedArchive: true,
  maxSize: "20m",
  maxFiles: "14d",
});

const accessTransport = new transports.DailyRotateFile({
  filename: config.accessLog,
  format: combine(
    timestamp(),
    printf((info) => `${info.level}: ${[info.timestamp]}: ${info.message}`)
  ),
  datePattern: "YYYY-MM-DD",
  zippedArchive: true,
  maxSize: "20m",
  maxFiles: "14d",
});

const logger = createLogger({
  transports: [errorTransport, accessTransport],
});

const logError = ({ req, message }) => {
  let url, method, httpVersion, headers;
  if (req) {
    ({ url, method, httpVersion, headers } = req);
  }
  logger.error({
    url,
    method,
    httpVersion,
    ...(headers && { host: headers.host, userAgent: headers["user-agent"] }),
    message,
  });
};

const logAccess = function ({ req, code, size }) {
  let url, method, httpVersion, headers;
  if (req) {
    ({ url, method, httpVersion, headers } = req);
  }

  logger.info({
    message: `${url} ${method} ${httpVersion} ${headers.host} ${headers["user-agent"]}`,
  });
};

module.exports = {
  logger,
  logAccess,
  logError,
};
