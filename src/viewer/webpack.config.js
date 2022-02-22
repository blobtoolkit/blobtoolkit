const path = require("path");
const webpack = require("webpack");
const PACKAGE = require("./package.json");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const { GitRevisionPlugin } = require("git-revision-webpack-plugin");
// const GitRevisionPlugin = require("git-revision-webpack-plugin");
const main = require("./src/config/main");
const CopyWebpackPlugin = require("copy-webpack-plugin");
const HtmlWebpackPlugin = require("html-webpack-plugin");

const devMode = process.env.NODE_ENV !== "production";

const BUILD_DIR = path.resolve(__dirname, "dist/public");
const APP_DIR = path.resolve(__dirname, "src/client/views");

const gitRevisionPlugin = new GitRevisionPlugin();

const protocol = main.https ? "https" : "http";

const config = {
  mode: devMode ? "development" : "production",
  entry: {
    main: ["@babel/polyfill", APP_DIR + "/index.jsx"],
  },
  output: {
    publicPath: main.mode == "production" ? main.basename + "/" : "/",
    path: BUILD_DIR + "/",
    filename: devMode ? "js/[name]/bundle.js" : "js/[name].[contenthash].js",
    chunkFilename: "js/[id].js",
  },
  resolve: {
    extensions: [".js", ".jsx"],
  },
  optimization: {
    splitChunks: {
      chunks: "all",
      maxInitialRequests: Infinity,
      minSize: 50000,
      cacheGroups: {
        defaultVendors: {
          test: /[\\/]node_modules[\\/]/,
          name(module) {
            const packageName = module.context.match(
              /[\\/]node_modules[\\/](.*?)([\\/]|$)/
            )[1];
            return `npm.${packageName.replace("@", "")}`;
          },
          chunks: "all",
        },
      },
    },
  },
  devServer: {
    hot: false,
    historyApiFallback: true,
    host: main.hostname,
    allowedHosts: "all",
    static: {
      directory: BUILD_DIR,
    },
    compress: true,
    port: main.client_port,
    // proxy: {
    //   "/api/**": { target: main.apiUrl },
    // },
  },
  devtool: "source-map",
  plugins: [
    new MiniCssExtractPlugin({
      filename: devMode ? "css/styles.css" : "css/[name].[contenthash].css",
      chunkFilename: "css/[id].css",
    }),
    gitRevisionPlugin,
    new webpack.DefinePlugin({
      API_URL: JSON.stringify(main.apiUrl),
      VERSION: JSON.stringify(PACKAGE.version),
      BASENAME: JSON.stringify(main.basename),
      STATIC_THRESHOLD: JSON.stringify(main.staticThreshold),
      NOHIT_THRESHOLD: JSON.stringify(main.nohitThreshold),
      CIRCLE_LIMIT: JSON.stringify(main.circleLimit),
      ABOUT: JSON.stringify(main.aboutUrl),
      HOME: JSON.stringify(protocol + "://" + main.hostname),
      GIT_VERSION: JSON.stringify(gitRevisionPlugin.version()),
      COMMIT_HASH: JSON.stringify(gitRevisionPlugin.commithash()),
      BRANCH: JSON.stringify(gitRevisionPlugin.branch()),
      GA_ID: JSON.stringify(main.ga_id),
      GDPR_URL: JSON.stringify(main.gdpr_url),
      DATASET_TABLE: main.dataset_table ? true : false,
      TARGET_TREE: main.target_tree ? true : false,
      MESSAGE: JSON.stringify(main.message),
    }),
    new HtmlWebpackPlugin({
      hash: true,
      title: "BlobToolKit - Viewer",
      template: "./src/client/index.html",
      minify: {
        collapseInlineTagWhitespace: true,
        collapseWhitespace: true,
        preserveLineBreaks: true,
        minifyURLs: true,
        removeComments: false,
        removeAttributeQuotes: true,
      },
    }),
    new CopyWebpackPlugin({
      patterns: [{ from: "./src/client/favicon" }],
    }),
    // new webpack.ExtendedAPIPlugin(),
  ],
  module: {
    rules: [
      {
        test: /\.(js|jsx)$/,
        include: APP_DIR,
        exclude: /(node_modules)/,
        use: {
          loader: "babel-loader",
          options: {
            presets: ["@babel/preset-env"],
          },
        },
      },
      {
        test: /\.html$/,
        include: APP_DIR,
        use: [
          {
            loader: "html-loader",
          },
        ],
      },
      {
        test: /\.css$/,
        use: [
          {
            loader: "style-loader",
            options: { injectType: "singletonStyleTag" },
          },
          {
            loader: "css-loader",
            options: {
              modules: false,
            },
          },
        ],
        include: [
          /node_modules/,
          path.resolve(__dirname, "src/client/views/style/node_modules.css"),
        ],
      },
      {
        test: /\.(sa|sc|c)ss$/,
        exclude: [
          /node_modules/,
          path.resolve(
            __dirname,
            "/src/client/views/components/style/node_modules.css"
          ),
        ],
        use: [
          MiniCssExtractPlugin.loader,
          {
            loader: "css-loader",
            options: {
              modules: {
                localIdentName: "[name]__[local]___[hash:base64:5]",
              },
              sourceMap: true,
              importLoaders: 2,
            },
          },
          "postcss-loader",
          "sass-loader",
        ],
      },
      {
        test: /\.svg$/,
        use: ["svg-sprite-loader", "svgo-loader"],
      },
      {
        test: /\.(gif|png|jpe?g)$/i,
        loader: "file-loader",
        options: {
          name: "img/[contenthash].[ext]",
          publicPath: main.mode == "production" ? main.basename + "/" : "/",
        },
      },
    ],
  },
};

module.exports = config;
