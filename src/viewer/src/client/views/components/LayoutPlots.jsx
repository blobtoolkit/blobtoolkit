import { getCatAxis, getMainPlot } from "../reducers/plot";
import {
  getDatasetID,
  getHashString,
  getQueryString,
  getStatic,
  getView,
  setQueryString,
  toggleHash,
} from "../reducers/location";
import {
  getDatasetIsActive,
  getNohitThreshold,
  getStaticThreshold,
} from "../reducers/repository";
import { getDatasetName, getSelectedDatasetMeta } from "../reducers/dataset";

import BuscoPlot from "./BuscoPlot";
import CumulativePlot from "./CumulativePlot";
import DatasetSpinner from "./DatasetSpinner";
import DatasetTable from "./DatasetTable";
import DetailPlot from "./DetailPlot";
import FindDatasets from "./FindDatasets";
import GetStarted from "./GetStarted";
import HomePage from "./HomePage";
import MainPlot from "./MainPlot";
import { NoHitWarning } from "./NoHitWarning";
import NoPlot from "./NoPlot";
import React from "react";
import SelectWarning from "./SelectWarning";
import SnailPlot from "./SnailPlot";
import StaticPlot from "./StaticPlot";
import TablePlot from "./TablePlot";
import TreeMapPlot from "./TreeMapPlot";
import { connect } from "react-redux";
import { getBinsForCat } from "../reducers/field";
import { getFields } from "../reducers/field";
import { getPlotShape } from "../reducers/plotParameters";
import { getWindowBinsForCat } from "../reducers/plotData";
import qs from "qs";
import { queryToStore } from "../querySync";
import styles from "./Layout.scss";

const dataset_table = DATASET_TABLE || false;

class PlotsLayoutComponent extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      datasetId: this.props.datasetId,
      static: this.props.static,
      keys: false,
      cat: false,
      warn: false,
    };
  }

  componentDidMount() {}
  componentDidUpdate() {
    if (!this.props.static && this.state.static) {
      if (
        this.props.meta &&
        this.props.meta.records > this.props.staticThreshold
      ) {
        this.setState({ static: this.props.static });
      }
    } else if (this.props.static && !this.state.static) {
      this.setState({ static: this.props.static });
    } else if (this.props.queryString) {
      if (
        this.state.keys &&
        !this.state.warn &&
        qs.parse(this.props.queryString)[this.state.cat] == this.state.keys
      ) {
        this.setState({ warn: true });
      } else if (
        this.state.warn &&
        qs.parse(this.props.queryString)[this.state.cat] != this.state.keys
      ) {
        this.setState({ warn: false });
      }
    } else if (this.state.warn) {
      this.setState({ warn: false });
    }
    if (!this.props.static) {
      if (
        this.props.meta &&
        this.props.meta.records > this.props.nohitThreshold &&
        this.props.bins &&
        this.props.cat &&
        !this.state.keys &&
        !this.state.cat
      ) {
        let index = this.props.bins.findIndex((x) => x.id == "no-hit");
        if (index > -1) {
          let keys = this.props.bins[index].keys;
          let cat = `${this.props.cat}--Keys`;
          let qstr = `${this.props.cat}--Keys=${keys.join(",")}`;
          this.setState({ keys: keys.join(","), cat });
          this.props.updateStore(qstr);
          this.setState({ warn: true });
        }
      } else if (
        this.props.meta &&
        this.props.meta.records <= this.props.nohitThreshold &&
        this.state.keys &&
        this.state.cat &&
        this.state.warn
      ) {
        let qstr = `${this.props.cat}--Keys=&nohitThreshold=${this.props.nohitThreshold}`;
        this.props.updateStore(qstr);
        this.setState({ keys: false, cat: false, warn: false });
      }
    }
  }
  render() {
    if (!this.props.datasetId) {
      return (
        <div className={styles.fill_parent}>
          <HomePage />
        </div>
      );
    } else if (!this.props.active || this.props.active == "loading") {
      return (
        <div className={styles.fill_parent}>
          <DatasetSpinner />
        </div>
      );
    }
    let warning = this.state.warn && (
      <NoHitWarning nohitThreshold={this.props.nohitThreshold} />
    );
    if (
      !this.props.datasetId ||
      (dataset_table && this.props.activeTab == "Datasets")
    ) {
      return <HomePage toggleHash={this.props.toggleHash} />;
    }
    let defaultPlot = <MainPlot {...this.props} />;
    let problem;
    if (
      this.props.active &&
      this.props.active != "loading" &&
      Object.keys(this.props.plot.axes).length < 4
    ) {
      if (!this.props.plot.axes.cat) {
        defaultPlot = <SnailPlot {...this.props} warning="noCat" />;
        problem = "noCat";
      } else {
        defaultPlot = <CumulativePlot {...this.props} warning="noBlob" />;
        problem = "noBlob";
      }
    }
    let view;

    if (this.props.static) {
      switch (this.props.view) {
        case "detail":
          view = <DetailPlot {...this.props} />;
          break;
        default:
          view = <StaticPlot {...this.props} />;
          break;
      }
    } else if (problem) {
      if (this.props.view == "blob") {
        view = defaultPlot;
      } else if (
        (problem == "noBlob" || problem == "noCat") &&
        this.props.view == "cumulative"
      ) {
        view = defaultPlot;
      }
    }
    if (!view) {
      switch (this.props.view || "blob") {
        case "busco":
          view = <BuscoPlot {...this.props} />;
          break;
        case "cumulative":
          console.log("here");
          view = <CumulativePlot {...this.props} />;
          break;
        case "detail":
          view = <DetailPlot {...this.props} />;
          break;
        case "notfound":
          view = <NoPlot {...this.props} />;
          break;
        case "snail":
          view = <SnailPlot {...this.props} />;
          break;
        case "table":
          view = <TablePlot {...this.props} />;
          break;
        case "treemap":
          view = <TreeMapPlot {...this.props} />;
          break;
        case "report":
          view = (
            <div className={styles.fill_parent}>
              <div className={styles.quarter}>
                <MainPlot {...this.props} />
              </div>
              <div className={styles.quarter}>
                <CumulativePlot {...this.props} />
              </div>
              <div className={styles.quarter}>
                <SnailPlot {...this.props} />
              </div>
              <div className={styles.quarter}>
                <BuscoPlot {...this.props} />
              </div>
              <div className={styles.quarter}>
                <DetailPlot {...this.props} />
              </div>
            </div>
          );
          break;
        default:
          view = defaultPlot;
          break;
      }
    }
    return (
      <div className={styles.fill_parent}>
        {view}
        <DatasetSpinner />
        {warning}
      </div>
    );
  }
}

class LayoutPlots extends React.Component {
  constructor(props) {
    super(props);
    (this.mapStateToProps = (state) => {
      let isStatic = getStatic(state);
      if (isStatic) {
        return {
          active: getDatasetIsActive(state),
          datasetId: getDatasetID(state),
          view: getView(state),
          plot: getMainPlot(state),
          activeTab: getHashString(state),
          static: true,
        };
      }
      let active = getDatasetIsActive(state);
      let cat, bins, view;
      if (active && active != "loading") {
        let fields = getFields(state);
        let shape = getPlotShape(state);
        cat = getCatAxis(state);
        view = getView(state);
        if (
          view == "blob" &&
          shape == "lines" &&
          fields &&
          fields[`${cat}_windows`]
        ) {
          bins = getWindowBinsForCat(state);
        } else {
          bins = getBinsForCat(state);
        }
      }
      return {
        active,
        datasetId: getDatasetID(state),
        datasetName: getDatasetName(state),
        plot: getMainPlot(state),
        activeTab: getHashString(state),
        view,
        bins,
        cat,
        staticThreshold: getStaticThreshold(state),
        nohitThreshold: getNohitThreshold(state),
        meta: getSelectedDatasetMeta(state),
        static: false,
        queryString: getQueryString(state),
      };
    }),
      (this.mapDispatchToProps = (dispatch) => {
        return {
          toggleHash: (value) => dispatch(toggleHash(value)),
          updateStore: (str, searchReplace) => {
            let values = qs.parse(str.replace("?", ""));
            dispatch(queryToStore({ values, searchReplace }));
          },
          updateQueryString: (qStr) => dispatch(setQueryString(qStr)),
        };
      });
  }

  render() {
    const ConnectedLayout = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(PlotsLayoutComponent);
    return <ConnectedLayout {...this.props} />;
  }
}

export default LayoutPlots;
