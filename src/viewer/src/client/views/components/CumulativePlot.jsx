import { fillParent, plotPaths } from "../reducers/plotStyles";
import { getCurveOrigin, getShowTotal } from "../reducers/plotParameters";

import AxisTitle from "./AxisTitle";
import CumulativePlotBoundary from "./CumulativePlotBoundary";
import { ExportButton } from "./ExportButton";
import FigureCaption from "./FigureCaption";
import { NoBlobWarning } from "./NoBlobWarning";
import PlotAxisTitle from "./PlotAxisTitle";
import PlotLegend from "./PlotLegend";
import React from "react";
import { connect } from "react-redux";
import { cumulativeCurves } from "../reducers/summary";
import { getSearchTerm } from "../reducers/location";
import { getSelectedDatasetMeta } from "../reducers/dataset";
import { getShowReference } from "../reducers/reference";
import styles from "./Plot.scss";

class Cumulative extends React.Component {
  render() {
    if (
      !this.props.cumulative ||
      !this.props.cumulative.zAxis ||
      this.props.cumulative.paths.byCat.length == 0
    )
      return null;
    let side = 1110;
    let viewbox = "0 0 " + side + " " + side;
    let legend = (
      <g transform="translate(685,705)">
        <PlotLegend />
      </g>
    );
    let colors = this.props.cumulative.palette.colors;
    let all = this.props.cumulative.paths.all;
    let yValues = this.props.cumulative.values.all;
    let records = this.props.cumulative.records;
    let span = this.props.cumulative.span;
    let byCat = this.props.cumulative.paths.byCat;
    let reference = this.props.cumulative.paths.reference;
    let searchTerm = this.props.searchTerm;
    let referenceIds = this.props.cumulative.referenceIds;
    let transform = "translate(0,0)";
    let yLabel = "cumulative " + this.props.cumulative.zAxis;
    let xLabel = "cumulative count";
    let paths = byCat.map((d, i) => {
      if (this.props.origin == "y") {
        let offsets = this.props.cumulative.paths.offsets;
        transform = "translate(" + offsets[i].x + "," + -offsets[i].y + ")";
      }
      if (this.props.origin == "x") {
        let offsets = this.props.cumulative.paths.count_offsets;
        transform = "translate(" + offsets[i].x + "," + -offsets[i].y + ")";
      }
      return (
        <path
          style={this.props.plotPaths.bold}
          d={d}
          key={i}
          fill="none"
          stroke={colors[i]}
          transform={transform}
          strokeLinecap="round"
        />
      );
    });
    let refPaths;
    let refText;
    if (this.props.showReference) {
      refPaths = reference.map((d, i) => {
        return (
          <path
            d={d}
            key={i}
            fill="none"
            stroke="#999"
            strokeWidth={1.5}
            strokeDasharray={"3 6"}
            strokeLinecap="round"
          />
        );
      });
      refText = referenceIds.map((obj, i) => {
        return (
          <text
            key={i}
            fill="#999"
            fontSize={16}
            style={{ cursor: "default", pointerEvents: "auto" }}
            textAnchor="end"
            transform={"translate(" + obj.offset.x + "," + obj.offset.y + ")"}
          >
            <a
              href={`/view/${searchTerm}/dataset/${obj.id}/cumulative#Settings`}
              target="blank"
              style={{ fill: "#999" }}
            >
              {obj.id}
            </a>
          </text>
        );
      });
    }

    let exportButtons = (
      <span className={styles.download}>
        <ExportButton
          view="cumulative"
          element="cumulative_plot"
          prefix={this.props.datasetId + ".cumulative"}
          format="svg"
        />
        <ExportButton
          view="cumulative"
          element="cumulative_plot"
          prefix={this.props.datasetId + ".cumulative"}
          format="png"
          size={side}
        />
      </span>
    );
    let warning;
    if (this.props.warning == "noBlob") {
      warning = <NoBlobWarning source="Cumulative" />;
    }
    return (
      <div className={styles.outer}>
        <div className={styles.fill_parent}>
          <svg
            id="cumulative_plot"
            ref={(elem) => {
              this.svg = elem;
            }}
            className={styles.main_plot + " " + styles.fill_parent}
            style={{ fontSize: "14px" }}
            viewBox={viewbox}
            preserveAspectRatio="xMinYMin"
          >
            <g transform={"translate(100,10)"}>
              <CumulativePlotBoundary
                yValues={yValues}
                records={records}
                span={span}
              />
              {refPaths}
              {refText}
              {this.props.showTotal && (
                <path
                  style={this.props.plotPaths.axis}
                  d={all}
                  fill="none"
                  stroke="#999"
                  strokeLinecap="round"
                />
              )}
              {paths}
              {legend}
              <AxisTitle axis="y" title={yLabel} side={side} />
              <AxisTitle axis="x" title={xLabel} side={side} />
            </g>
          </svg>
          {exportButtons}
        </div>
        {warning}
        <FigureCaption {...this.props} />
      </div>
    );
  }
}

class CumulativePlot extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      return {
        cumulative: cumulativeCurves(state),
        meta: getSelectedDatasetMeta(state, this.props.datasetId),
        origin: getCurveOrigin(state),
        showTotal: getShowTotal(state),
        searchTerm: getSearchTerm(state),
        showReference: getShowReference(state),
        plotPaths: plotPaths(state),
      };
    };
  }

  render() {
    const ConnectedCumulative = connect(this.mapStateToProps)(Cumulative);
    return <ConnectedCumulative {...this.props} />;
  }
}

export default CumulativePlot;
