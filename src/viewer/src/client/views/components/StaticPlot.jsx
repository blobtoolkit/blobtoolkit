import { getPlotShape, getPngResolution } from "../reducers/plotParameters";

import { ExportButton } from "./ExportButton";
import { NoBlobWarning } from "./NoBlobWarning";
import React from "react";
import StaticBuscoPlot from "./StaticBuscoPlot";
import StaticMissing from "./StaticMissing";
import StaticWarning from "./StaticWarning";
import { apiUrl } from "../reducers/api";
import { connect } from "react-redux";
import { getSelectedDatasetMeta } from "../reducers/dataset";
import { getStaticFields } from "../reducers/dataset";
import { getStaticThreshold } from "../reducers/repository";
import styles from "./Plot.scss";

const arrayBufferToBase64 = (buffer) => {
  let binary = "";
  let bytes = [].slice.call(new Uint8Array(buffer));
  bytes.forEach((b) => (binary += String.fromCharCode(b)));
  return window.btoa(binary);
};

class Static extends React.Component {
  constructor(props) {
    super(props);
    this.state = { source: null, available: true, status: 0 };
  }

  getResolution() {
    let container = this.refs.static_image.parentNode;
    let res = Math.min(container.clientWidth, container.clientHeight);
    let query =
      "(-webkit-min-device-pixel-ratio: 2), (min-device-pixel-ratio: 2), (min-resolution: 192dpi)";
    if (window.matchMedia(query).matches) {
      res *= 2;
    }
    return res;
  }

  componentDidMount() {
    if (this.refs.static_image && this.props.hasStatic) {
      let shape;
      let view = this.props.view;
      if (view == "blob" && Object.keys(this.props.plot.axes).length < 4) {
        view = "cumulative";
      } else if (view == "blob") {
        shape = this.props.shape;
      }
      this.fetchImage(view, shape, this.getResolution());
    }
  }
  componentDidUpdate(nextProps) {
    if (
      this.refs.static_image &&
      this.props.hasStatic &&
      nextProps.shape != this.props.shape
    ) {
      let shape;
      if (this.props.view == "blob") {
        shape = nextProps.shape;
      }
      this.fetchImage(this.props.view, shape, this.getResolution());
    }
  }

  imageURL(id, view, shape, width, format = "png") {
    let url = apiUrl + "/image/" + id + "/" + view;
    if (shape) {
      url += "/" + shape;
    }
    url += "?format=" + format;
    if (width) {
      url += "&width=" + width;
    }
    return url;
  }

  fetchImage(view, shape, width) {
    width = width || this.props.pngResolution;
    let url;
    if (this.state.status == 0) {
      url = this.imageURL(this.props.datasetId, view, shape, width / 4, "png");
    }
    if (this.state.status == 1) {
      url = this.imageURL(this.props.datasetId, view, shape, width, "png");
    }
    fetch(url).then(
      (response) => {
        if (!response.ok) {
          this.setState({ available: false });
          return;
        }
        response.arrayBuffer().then((buffer) => {
          let base64Flag = "data:image/png;base64,";
          let imageStr = arrayBufferToBase64(buffer);
          if (this.refs.static_image) {
            this.refs.static_image.src = base64Flag + imageStr;
            this.setState({ available: true });
            if (this.state.status == 0) {
              this.setState({ status: 1 });
              this.fetchImage(view, shape, width);
            } else {
              this.setState({ image: buffer });
            }
          }
        });
      },
      (error) => console.log("An error occured.", error)
    );
  }

  render() {
    let view = this.props.view;
    let warning;
    if (view == "blob" && Object.keys(this.props.plot.axes).length < 4) {
      warning = <NoBlobWarning source="Cumulative" />;
      view = "cumulative";
    } else if (this.state.available) {
      warning = (
        <StaticWarning
          name={this.props.meta.name}
          threshold={this.props.threshold}
          records={this.props.meta.records}
        />
      );
    } else {
      warning = (
        <StaticMissing name={this.props.meta.name} view={this.props.view} />
      );
    }
    let shape;
    let prefix = this.props.datasetId + "." + view;
    if (view == "blob") {
      shape = this.props.shape;
      prefix += "." + shape;
    }
    let exportButtons = (
      <span className={styles.download}>
        <ExportButton
          url={this.imageURL(
            this.props.datasetId,
            view,
            shape,
            this.props.pngResolution,
            "svg"
          )}
          prefix={prefix}
          format="svg"
        />
        <ExportButton
          url={this.imageURL(
            this.props.datasetId,
            view,
            shape,
            this.props.pngResolution,
            "png"
          )}
          prefix={prefix}
          format="png"
        />
      </span>
    );
    if (view == "busco") {
      return (
        <div className={styles.fill_parent}>
          <StaticBuscoPlot />
          {warning}
        </div>
      );
    }
    return (
      <div className={styles.fill_parent + " " + styles.centered_content}>
        <div className={styles.outer}>
          <img
            id="static_image"
            className={styles.static_image}
            ref="static_image"
          />
          {exportButtons}
        </div>
        {warning}
      </div>
    );
  }
}

class StaticPlot extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      return {
        hasStatic: getStaticFields(state),
        meta: getSelectedDatasetMeta(state),
        threshold: getStaticThreshold(state),
        shape: getPlotShape(state),
        pngResolution: getPngResolution(state),
      };
    };
    this.mapDispatchToProps = (dispatch) => {
      return {};
    };
  }

  render() {
    const ConnectedStatic = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(Static);
    return <ConnectedStatic {...this.props} />;
  }
}

export default StaticPlot;
