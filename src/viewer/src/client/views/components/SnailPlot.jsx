import { addRecords, removeRecords, setSelectSource } from "../reducers/select";
import {
  chooseCircumferenceScale,
  chooseRadiusScale,
  getCurveOrigin,
  getLargeFonts,
} from "../reducers/plotParameters";
import {
  circularCurves,
  circularSelection,
  getBuscoData,
  getBuscoPaths,
  getBuscoSets,
  getCircular,
} from "../reducers/summary";
import { fetchRawData, fields, getDetailsForFieldId } from "../reducers/field";
import { fillParent, plotPaths, plotText } from "../reducers/plotStyles";

import { CircleAxis } from "./CircleAxis";
import { ExportButton } from "./ExportButton";
import FigureCaption from "./FigureCaption";
import Pointable from "react-pointable";
import React from "react";
import SnailPlotLegend from "./SnailPlotLegend";
import SnailPlotScale from "./SnailPlotScale";
import { connect } from "react-redux";
import { cumulativeCurves } from "../reducers/summary";
import { arc as d3Arc } from "d3-shape";
import { format as d3format } from "d3-format";
import { scaleLinear as d3scaleLinear } from "d3-scale";
import { getColorScheme } from "../reducers/color";
import { getSelectedDatasetMeta } from "../reducers/dataset";
import styles from "./Plot.scss";

const xScale = d3scaleLinear().range([-425, 425]);
const yScale = d3scaleLinear().range([425, -425]);
const radialCoords = (event) => {
  let bounds = event.target.getBoundingClientRect();
  xScale.domain([bounds.left, bounds.right]);
  yScale.domain([bounds.top, bounds.bottom]);
  let x = xScale(event.clientX);
  let y = yScale(event.clientY);
  let theta = Math.atan2(x, y);
  theta = theta >= 0 ? theta : 2 * Math.PI - Math.abs(theta);
  let h = Math.sqrt(x ** 2 + y ** 2);
  return { theta, h };
};

const SegmentStats = ({
  theta,
  h,
  ratio,
  overlay,
  sums,
  nXlens,
  nXnums,
  gcs,
  ns,
  index,
}) => {
  if (!overlay) return {};
  let cScale = d3scaleLinear()
    .domain([0, 2 * Math.PI * ratio])
    .range([0, 100]);
  let start = Math.floor(cScale(theta));
  let end = start + 1;
  let gc = {},
    allGC = [],
    allN = [],
    nX = {},
    path,
    N = {},
    indices = [];
  if (start < 100) {
    let innerRadius = 0;
    let outerRadius = 350;
    if (h >= 350) {
      innerRadius = 350;
      outerRadius = 425;
      for (let i = start * 10; i < end * 10; i++) {
        allGC = allGC.concat(gcs[i].gc);
        allN = allN.concat(ns[i].n);
        indices = indices.concat(index[i]);
      }
      N.min = Math.min(...allN);
      N.max = Math.max(...allN);
      N.count = allN.length;
      N.mean = allN.reduce((a = 0, b) => a + b) / N.count;
      gc.min = Math.min(...allGC);
      gc.max = Math.max(...allGC);
      gc.count = allGC.length;
      gc.mean = allGC.reduce((a = 0, b) => a + b) / gc.count;
    } else {
      start = 0;
      nX.len = nXlens[end * 10 - 1];
      nX.num = nXnums[end * 10 - 1];
    }
    path = d3Arc()({
      startAngle: (Math.PI * ratio * start) / 50,
      endAngle: (Math.PI * ratio * end) / 50,
      innerRadius,
      outerRadius,
    });
  }
  return { pct: end, path, gc, N, nX, indices };
};

const SnailSegment = ({ path }) => (
  <path d={path} fill={"rgba(225,225,225,0.4)"} stroke={"none"} />
);

class Snail extends React.Component {
  constructor(props) {
    super(props);
    // let ratio = 1
    if (!props.data) {
      return null;
    }
    let ratio = props.selection.ratio;

    this.state = {
      theta: 3.14,
      h: 320,
      ratio,
      overlay: false,
      mouseDown: false,
      last: -1,
      action: false,
    };
  }

  toggleSelect(theta, h) {
    let ratio = this.state.ratio;
    let cScale = d3scaleLinear()
      .domain([0, 2 * Math.PI * ratio])
      .range([0, 100]);
    let cell = Math.floor(cScale(theta));
    let action = this.state.action;
    if (cell >= 0 && cell < 100 && cell != this.state.last) {
      if (h <= 450) {
        let indices = [];
        let i = h < 350 ? 0 : cell * 10;
        for (i; i < (cell + 1) * 10; i++) {
          indices = indices.concat(this.props.data.values.index[i]);
          indices = [...new Set(indices)];
        }
        if (!action) {
          action = this.props.selection.segments[cell].selected
            ? "remove"
            : "add";
        }
        if (action == "remove") {
          this.props.removeRecords(indices);
          this.props.setSelectSource();
        } else if (action == "add") {
          this.props.addRecords(indices);
          this.props.setSelectSource();
        }
      }
      this.setState({ last: cell, action });
    }
  }

  componentDidMount() {
    let fields = ["gc", "length"];
    if (this.props.ncountMeta.meta.datatype) {
      fields.push("ncount");
    }
    fields.forEach((f) => {
      if (!this.props[f]) {
        this.props.activate(f);
      }
    });
    if (this.props.buscoSets && !this.props.buscoData) {
      this.props.activate(this.props.buscoSets[0]);
    }
  }

  render() {
    if (!this.props.circular) return null;
    let format = d3format(".2s");
    let pctFormat = d3format(".1%");
    let commaFormat = d3format(",");
    let side = 1050;
    let viewbox = "0 0 " + side + " " + side;
    let pathProps = this.props.circular.pathProps;
    let paths = [];
    let legend = this.props.circular.legend;
    let bottomLeft = (
      <SnailPlotScale
        title={"Scale"}
        scale={this.props.circular.scale}
        onChangeCircumference={this.props.onChangeCircumference}
        onChangeRadius={this.props.onChangeRadius}
      />
    );
    let segment;
    let sums = this.props.data.values.sum;
    let nXlens = this.props.data.values.nXlen;
    let nXnums = this.props.data.values.nXnum;
    let gcs = this.props.data.values.gc;
    let ns = this.props.data.values.n;
    let index = this.props.data.values.index;
    let largeFonts = this.props.largeFonts;
    let stats = SegmentStats({
      ...this.state,
      sums,
      gcs,
      ns,
      nXlens,
      nXnums,
      index,
    });
    let topLeft, bottomRight, topRight;
    if (legend.composition) {
      bottomRight = (
        <SnailPlotLegend title={"Composition"} list={legend.composition} />
      );
    }
    if (stats.path) {
      segment = <SnailSegment path={stats.path} />;
      let list = legend.composition.map((l) => ({
        label: l.label,
        value: l.value,
        color: l.color,
      }));
      if (stats.gc && stats.gc.mean) {
        let gcs = [stats.gc.mean, stats.gc.min, stats.gc.max];
        let ns = [stats.N.mean, stats.N.min, stats.N.max];
        list[0].label = pctFormat(gcs[0]);
        list[0].value = pctFormat(gcs[1]) + "–" + pctFormat(gcs[2]);
        list[1].label = pctFormat(1 - gcs[0]);
        list[1].value = pctFormat(1 - gcs[2]) + "–" + pctFormat(1 - gcs[1]);
        if (list[2]) {
          list[2].label = pctFormat(ns[0]);
          list[2].value = pctFormat(ns[1]) + "–" + pctFormat(ns[2]);
        }
      }
      bottomRight = <SnailPlotLegend title={"Composition"} list={list} />;
      if (stats.nX && stats.nX.len) {
        let list = [
          {
            label: "count: " + commaFormat(stats.nX.num),
            color: legend.stats[0].color,
          },
          {
            label: "length: " + commaFormat(stats.nX.len),
            color: legend.stats[1].color,
          },
        ];
        bottomRight = <SnailPlotLegend title={"N" + stats.pct} list={list} />;
      }
    }
    if (legend.stats) {
      topLeft = (
        <SnailPlotLegend title={"Scaffold statistics"} list={legend.stats} />
      );
    }
    if (this.props.buscoPaths) {
      let buscoListA = [];
      let buscoListB = [];
      let subtitle = [];
      buscoListA.push({
        label: largeFonts ? "Comp." : "Complete",
        value: pctFormat(this.props.buscoData.fractions.c),
        color: "rgb(51, 160, 44)",
      });
      buscoListA.push({
        label: largeFonts ? "Dupl." : "Duplicated",
        value: pctFormat(this.props.buscoData.fractions.d),
        color: "rgb(32, 100, 27)",
      });
      if (this.props.buscoMeta.meta.set) {
        subtitle.push({
          label: this.props.buscoMeta.meta.set,
          value: this.props.buscoData.total,
        });
      } else {
        subtitle.push({ label: "Total: " + this.props.buscoData.total });
      }
      buscoListB.push({
        label: largeFonts ? "Frag." : "Fragmented",
        value: pctFormat(this.props.buscoData.fractions.f),
        color: "rgb(163, 226, 127)",
      });
      buscoListB.push({
        label: "Missing",
        value: pctFormat(this.props.buscoData.fractions.m),
        color: "white",
      });

      let axis = this.props.circular.axes.inner;
      let buscoAxis = <CircleAxis radius={200} />;

      topRight = (
        <g transform="translate(10,0)">
          <g transform="translate(-225,0)">
            <SnailPlotLegend title={"BUSCO"} list={buscoListA} />
          </g>
          <g transform="translate(-110,-25)">
            <SnailPlotLegend title={""} list={subtitle} />
          </g>
          <g transform="translate(-20,0)">
            <SnailPlotLegend title={""} list={buscoListB} />
          </g>
          <g transform="translate(70,170) scale(0.35)">
            <path
              d={this.props.buscoPaths.c}
              fill={"rgb(51, 160, 44)"}
              stroke={"none"}
            />
            <path
              d={this.props.buscoPaths.d}
              fill={"rgb(32, 100, 27)"}
              stroke={"none"}
            />
            <path
              d={this.props.buscoPaths.f}
              fill={"rgb(163, 226, 127)"}
              stroke={"none"}
            />
            {buscoAxis}
          </g>
        </g>
      );
    }
    Object.keys(this.props.circular.paths).forEach((k, i) => {
      let d = this.props.circular.paths[k];
      paths.push(
        <path
          d={d}
          key={k}
          fill={pathProps[k].fill}
          stroke={pathProps[k].stroke}
          strokeWidth={pathProps[k].strokeWidth || 2}
          strokeDasharray={pathProps[k].strokeDasharray || null}
          strokeLinecap="round"
        />
      );
    });
    let axes = this.props.circular.axes;
    Object.keys(axes).forEach((k, i) => {
      let axis = axes[k];
      paths.push(
        <path
          style={this.props.plotPaths.axis}
          d={axis.path}
          key={k}
          fill="none"
          stroke="black"
          strokeLinecap="round"
        />
      );
      if (axis.ticks) {
        axis.ticks.major.forEach((d, idx) => {
          paths.push(
            <path
              style={this.props.plotPaths.axis}
              d={d}
              key={k + "_major_" + idx}
              fill="none"
              stroke="black"
              strokeLinecap="round"
            />
          );
        });
        axis.ticks.minor.forEach((d, idx) => {
          paths.push(
            <path
              style={this.props.plotPaths.fine}
              d={d}
              key={k + "_minor_" + idx}
              fill="none"
              stroke="black"
              strokeLinecap="round"
            />
          );
        });
        if (axis.ticks.labels) {
          axis.ticks.labels.forEach((d, idx) => {
            let textAnchor = "middle";
            let startOffset = "50%";
            let fontSize = d.fontSize ? d.fontSize : "18px";
            if (d.align == "left") {
              textAnchor = "left";
              startOffset = 0;
            }
            if (d.align == "right") {
              textAnchor = "end";
              startOffset = "100%";
            }
            paths.push(
              <path
                style={this.props.plotPaths.axis}
                d={d.path}
                key={k + "_path_" + idx}
                id={k + "_path_" + idx}
                fill="none"
                stroke="none"
              />
            );
            paths.push(
              <text
                key={k + "_text_" + idx}
                style={Object.assign({}, this.props.plotText.axisLabel, {
                  fontSize,
                })}
              >
                <textPath
                  xlinkHref={"#" + k + "_path_" + idx}
                  textAnchor={textAnchor}
                  startOffset={startOffset}
                >
                  {d.text}
                </textPath>
              </text>
            );
          });
        }
      }
    });
    this.props.selection.paths.forEach((path, i) => {
      paths.push(
        <path
          d={path.path}
          key={`selection_${i}`}
          fill={path.partial ? "none" : this.props.colors.halfHighlightColor}
          stroke={this.props.colors.highlightColor}
          strokeWidth={3}
          strokeLinecap="round"
        />
      );
    });

    let exportButtons = (
      <span className={styles.download}>
        <ExportButton
          view="snail"
          element="snail_plot"
          prefix={this.props.datasetId + ".snail"}
          format="svg"
        />
        <ExportButton
          view="snail"
          element="snail_plot"
          prefix={this.props.datasetId + ".snail"}
          format="png"
          size={side}
        />
      </span>
    );
    return (
      <div className={styles.outer}>
        <div className={styles.fill_parent}>
          <svg
            id="snail_plot"
            ref={(elem) => {
              this.svg = elem;
            }}
            className={styles.main_plot + " " + styles.fill_parent}
            style={{ fontSize: "14px" }}
            viewBox={viewbox}
            preserveAspectRatio="xMinYMin"
          >
            <g transform={"translate(10,1000)"}>
              <text style={this.props.plotText.snailPlotTitle}>
                Dataset: {this.props.meta.id}
              </text>
            </g>
            <g transform={"translate(525,550)"}>
              {paths}
              {segment}
            </g>
            <Pointable
              tagName="g"
              onPointerMove={(e) => {
                e.preventDefault();
                let coords = radialCoords(e);
                if (this.state.mouseDown) {
                  this.toggleSelect(coords.theta, coords.h);
                }
                this.setState({ ...coords, overlay: true });
              }}
              onPointerLeave={(e) => {
                e.preventDefault();
                this.setState({
                  overlay: false,
                  mouseDown: false,
                  last: -1,
                  action: false,
                });
              }}
              onPointerDown={(e) => {
                e.preventDefault();
                let coords = radialCoords(e);
                this.toggleSelect(coords.theta, coords.h);
                this.setState({ ...coords, overlay: true, mouseDown: true });
              }}
              onPointerUp={(e) => {
                e.preventDefault();
                this.setState({
                  overlay: false,
                  mouseDown: false,
                  last: -1,
                  action: false,
                });
              }}
            >
              <circle
                r={450}
                cx={520}
                cy={550}
                fill="rgba(255,255,255,0)"
                stroke="none"
                style={{
                  pointerEvents: "auto",
                  stroke: "none",
                  cursor: "pointer",
                }}
              />
            </Pointable>
            <g transform={"translate(10,35)"}>{topLeft}</g>
            <g transform={"translate(10,890)"}>{bottomLeft}</g>
            <g transform={"translate(850,890)"}>{bottomRight}</g>
            <g transform={"translate(850,35)"}>{topRight}</g>
          </svg>
          {exportButtons}
        </div>
        <FigureCaption {...this.props} />
      </div>
    );
  }
}

class SnailPlot extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      let buscoSets = getBuscoSets(state);
      let buscoData, buscoPaths, buscoMeta;
      if (buscoSets && buscoSets.length > 0) {
        (buscoData = getBuscoData(state, buscoSets[0])),
          (buscoPaths = getBuscoPaths(state, buscoSets[0])),
          (buscoMeta = getDetailsForFieldId(state, buscoSets[0]));
      }
      let ncountMeta = getDetailsForFieldId(state, "ncount");
      return {
        data: getCircular(state),
        // gc: getRawDataForFieldId(state, "gc"),
        // length: getRawDataForFieldId(state, "length"),
        // ncount: getRawDataForFieldId(state, "ncount"),
        circular: circularCurves(state),
        meta: getSelectedDatasetMeta(state),
        selection: circularSelection(state),
        colors: getColorScheme(state),
        plotPaths: plotPaths(state),
        plotText: plotText(state),
        largeFonts: getLargeFonts(state),
        buscoSets,
        buscoData,
        buscoPaths,
        buscoMeta,
        ncountMeta,
      };
    };
    this.mapDispatchToProps = (dispatch) => {
      return {
        onChangeCircumference: (value) =>
          dispatch(chooseCircumferenceScale(value)),
        onChangeRadius: (value) => dispatch(chooseRadiusScale(value)),
        activate: (id) => {
          dispatch(fetchRawData(id));
        },
        addRecords: (arr) => dispatch(addRecords(arr)),
        removeRecords: (arr) => dispatch(removeRecords(arr)),
        setSelectSource: () => dispatch(setSelectSource("snail")),
      };
    };
  }

  render() {
    const ConnectedSnail = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(Snail);
    return <ConnectedSnail {...this.props} />;
  }
}

export default SnailPlot;
