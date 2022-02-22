import {
  addRecords,
  changeSelectPolygon,
  getSelectPolygon,
  getSelectSource,
  getSelectedRecords,
  getSelectionDisplay,
  replaceRecords,
  selectNone,
  setSelectSource,
} from "../reducers/select";
import {
  curveBasis as d3CurveBasis,
  curveBasisClosed as d3CurveBasisClosed,
  curveCardinal as d3CurveCardinal,
  curveCardinalClosed as d3CurveCardinalClosed,
  curveLinear as d3CurveLinear,
  curveLinearClosed as d3CurveLinearClosed,
  curveMonotoneX as d3CurveMonotoneX,
  curveMonotoneY as d3CurveMonotoneY,
  line as d3Line,
} from "d3-shape";
import { getCatAxis, getXAxis, getYAxis, getZAxis } from "../reducers/plot";
import {
  getErrorBars,
  getPlotScale,
  getWindowSize,
  getZScale,
} from "../reducers/plotParameters";

import React from "react";
import { connect } from "react-redux";
import { polygonHull as d3PolygonHull } from "d3-polygon";
import { fetchRawData } from "../reducers/field";
import { getLinesPlotData } from "../reducers/plotData";
import { plotShapes } from "../reducers/plotStyles";

const smoothLineX = d3Line()
  .x((d) => d[0])
  .y((d) => d[1])
  .curve(d3CurveMonotoneX);

const smoothLineY = d3Line()
  .x((d) => d[0])
  .y((d) => d[1])
  .curve(d3CurveMonotoneY);

const cardinalLine = d3Line()
  .x((d) => d[0])
  .y((d) => d[1])
  .curve(d3CurveCardinal);

const closedLine = d3Line()
  .x((d) => d[0])
  .y((d) => d[1])
  .curve(d3CurveBasisClosed);

export default class PlotLinesSVG extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => ({
      xAxis: getXAxis(state),
      yAxis: getYAxis(state),
      zAxis: getZAxis(state),
      catAxis: getCatAxis(state),
      zScale: getZScale(state),
      plotScale: getPlotScale(state),
      windowSize: getWindowSize(state),
      errorBars: getErrorBars(state),
      data: getLinesPlotData(state),
      showSelection: getSelectionDisplay(state),
      selectedRecords: getSelectedRecords(state),
    });
    this.mapDispatchToProps = (dispatch) => ({
      fetchData: (id) => dispatch(fetchRawData(id)),
      addRecords: (arr) => {
        dispatch(setSelectSource("lines"));
        dispatch(addRecords(arr));
      },
      replaceRecords: (arr) => {
        dispatch(setSelectSource("lines"));
        dispatch(replaceRecords(arr));
      },
      selectNone: () => {
        dispatch(setSelectSource("lines"));
        dispatch(selectNone());
      },
    });
  }

  render() {
    const ConnectedLines = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(LinesSVG);
    return <ConnectedLines {...this.props} />;
  }
}

class LinesSVG extends React.Component {
  componentDidMount() {
    if (
      !this.props.data ||
      !this.props.data.coords ||
      this.props.data.coords.length == 0
    ) {
      let windowSuffix =
        this.props.windowSize == 0.1 ? "" : `_${this.props.windowSize}`;
      this.props.fetchData(`${this.props.xAxis}_windows${windowSuffix}`);
      this.props.fetchData(`${this.props.yAxis}_windows${windowSuffix}`);
      this.props.fetchData(`${this.props.zAxis}_windows${windowSuffix}`);
      this.props.fetchData(`${this.props.catAxis}_windows${windowSuffix}`);
    }
  }

  // componentDidUpdate() {
  //   if (
  //     !this.props.data ||
  //     !this.props.data.coords ||
  //     this.props.data.coords.length == 0
  //   ) {
  //     this.props.fetchData(`${this.props.xAxis}_windows_1000000`);
  //     this.props.fetchData(`${this.props.yAxis}_windows_1000000`);
  //     this.props.fetchData(`${this.props.zAxis}_windows_1000000`);
  //     this.props.fetchData(`${this.props.catAxis}_windows_1000000`);
  //   }
  // }

  handleClick(e, id) {
    e.stopPropagation();
    e.preventDefault();
    let arr = [...this.props.selectedRecords];
    if (arr.includes(id)) {
      this.props.replaceRecords(arr.filter((val) => val != id));
    } else {
      this.props.addRecords([id]);
    }
  }

  render() {
    let coords = this.props.data.coords;
    if (!coords || coords.length == 0) {
      return null;
    }
    let colors = this.props.data.colors;
    let paths = [];
    let selection = [];
    let selectedById = {};
    let highlightColor = "rgb(156, 82, 139)";
    this.props.selectedRecords.forEach((id) => {
      selectedById[id] = true;
    });
    coords.forEach((group, i) => {
      if (group.x.length > 0) {
        let points = [];
        let groupCircles = [];
        let selCircles = [];
        group.x.forEach((x, j) => {
          points.push([x, group.y[j]]);
          let color;
          let opacity;
          let sel;
          if (isNaN(group.cats[j])) {
            color = "white";
            opacity = 0.3;
            sel = false;
          } else {
            color = colors[group.cats[j]] || colors[9];
            opacity = 0.6;
            sel = true;
          }
          if (sel && selectedById[group.id]) {
            selCircles.push(
              <circle
                key={`${group.id}_${j}_sel`}
                cx={x}
                cy={group.y[j]}
                r={group.r[j]}
                fill={color}
                fillOpacity={1}
                style={{
                  strokeWidth: this.props.showSelection ? "6px" : "2px",
                  stroke: this.props.showSelection
                    ? highlightColor
                    : "rgb(89, 101, 111)",
                }}
                onPointerDown={(e) => this.handleClick(e, group.id)}
              />
            );
          }
          groupCircles.push(
            <circle
              key={`${group.id}_${j}`}
              cx={x}
              cy={group.y[j]}
              r={group.r[j]}
              fill={color}
              fillOpacity={selectedById[group.id] ? 1 : opacity}
              style={{
                strokeWidth: "1px",
                stroke: "rgb(89, 101, 111)",
              }}
              onPointerDown={(e) => this.handleClick(e, group.id)}
            />
          );
        });
        groupCircles = groupCircles.concat(selCircles);

        let groupPaths = [];
        let mainPath = cardinalLine(points);

        if (selectedById[group.id]) {
          if (
            (this.props.errorBars && group.ysd && !group.xsd) ||
            (group.xsd && !group.ysd)
          ) {
            let upperPoints = [];
            let lowerPoints = [];
            let upperPath, lowerPath;

            if (group.ysd) {
              group.x.forEach((x, j) => {
                upperPoints.push([x, group.y[j] + group.ysd[j]]);
                lowerPoints.push([x, group.y[j] - group.ysd[j]]);
              });
              upperPath = smoothLineX(upperPoints);
              lowerPath = smoothLineX(lowerPoints.reverse());
            } else {
              group.x.forEach((x, j) => {
                upperPoints.push([x + group.xsd[j], group.y[j]]);
                lowerPoints.push([x - group.xsd[j], group.y[j]]);
              });
              upperPath = smoothLineY(upperPoints);
              lowerPath = smoothLineY(lowerPoints.reverse());
            }
            groupPaths.push(
              <path
                key={`${group.id}_bounds`}
                stroke={"none"}
                fill="black"
                fillOpacity={0.1}
                d={upperPath + lowerPath.replace(/^M/, "L")}
                onPointerDown={(e) => this.handleClick(e, group.id)}
              />
            );
            groupPaths.push(
              <path
                key={`${group.id}_upper`}
                style={{ strokeWidth: "5px" }}
                stroke={"black"}
                fill="none"
                strokeOpacity={0.4}
                d={upperPath}
                onPointerDown={(e) => this.handleClick(e, group.id)}
              />
            );
            groupPaths.push(
              <path
                key={`${group.id}_lower`}
                style={{ strokeWidth: "5px" }}
                stroke={"black"}
                fill="none"
                strokeOpacity={0.4}
                d={lowerPath}
                onPointerDown={(e) => this.handleClick(e, group.id)}
              />
            );
            // groupPaths.push(
            //   <polyline
            //     key={`${group.id}_lower`}
            //     style={{ strokeWidth: "10px" }}
            //     stroke={"black"}
            //     strokeLinejoin="round"
            //     fill="none"
            //     points={lowerPoints.join(" ")}
            //     onPointerDown={(e) => this.handleClick(e, group.id)}
            //   />
            // );
            // mainPath = midPath;
          } else if (group.xsd && group.ysd) {
            let sdPoints = [];
            group.x.forEach((x, j) => {
              sdPoints.push([x - group.xsd[j], group.y[j] + group.ysd[j]]);
              sdPoints.push([x + group.xsd[j], group.y[j] + group.ysd[j]]);
              sdPoints.push([x + group.xsd[j], group.y[j] - group.ysd[j]]);
              sdPoints.push([x - group.xsd[j], group.y[j] - group.ysd[j]]);
            });

            let boundary = closedLine(d3PolygonHull(sdPoints));

            groupPaths.push(
              <path
                key={`${group.id}_area`}
                stroke={"none"}
                fill="black"
                fillOpacity={0.1}
                d={boundary}
                onPointerDown={(e) => this.handleClick(e, group.id)}
              />
            );
            groupPaths.push(
              <path
                key={`${group.id}_hull`}
                style={{ strokeWidth: "5px" }}
                stroke={"black"}
                fill="none"
                strokeOpacity={0.4}
                strokeLinejoin="round"
                fill="none"
                d={boundary}
                onPointerDown={(e) => this.handleClick(e, group.id)}
              />
            );
          }

          {
            if (this.props.showSelection) {
              groupPaths.push(
                <path
                  key={`${group.id}_sel`}
                  style={{ strokeWidth: "14px" }}
                  stroke={highlightColor}
                  strokeLinejoin="round"
                  fill="none"
                  d={mainPath}
                  onPointerDown={(e) => this.handleClick(e, group.id)}
                />
              );
            } else {
              groupPaths.push(
                <path
                  key={`${group.id}_sel`}
                  style={{ strokeWidth: "8px" }}
                  stroke={"rgb(89, 101, 111)"}
                  strokeLinejoin="round"
                  fill="none"
                  d={mainPath}
                  onPointerDown={(e) => this.handleClick(e, group.id)}
                />
              );
            }
          }
        }
        groupPaths.push(
          <path
            key={group.id}
            style={{ strokeWidth: selectedById[group.id] ? "4px" : "2px" }}
            stroke={colors[group.cat]}
            strokeLinejoin="round"
            fill="none"
            d={mainPath}
            onPointerDown={(e) => this.handleClick(e, group.id)}
          />
        );
        if (selectedById[group.id]) {
          selection.push(...groupPaths, ...groupCircles);
        } else {
          paths.push(...groupPaths, ...groupCircles);
        }
      }
      //   let kite = <g key={i}
      //                 style={{strokeWidth:"1px"}}
      //                 transform={`rotate(${coords[i].angle},${coords[i].y[0][0]},${coords[i].x[0][1]})`}
      //                 stroke={colors[i]}
      //                 fill="none">
      //                 <line key={`${i}_x`}
      //                       style={{strokeWidth:"3px"}}
      //                       x1={coords[i].x[0][0]}
      //                       y1={coords[i].x[0][1]}
      //                       x2={coords[i].x[1][0]}
      //                       y2={coords[i].x[1][1]}/>
      //                 <line key={`${i}_y`}
      //                       x1={coords[i].y[0][0]}
      //                       y1={coords[i].y[0][1]}
      //                       x2={coords[i].y[1][0]}
      //                       y2={coords[i].y[1][1]}/>
      //                 <polygon key={`${i}_poly`}
      //                          style={{strokeWidth:`3px`}}
      //                          points={coords[i].poly.map(c=>c[0]+','+c[1]).join(' ')}/>
      //               </g>
      //   paths.push( kite )
      // }
    });
    return (
      <g
        transform="translate(0, 0)"
        style={{ cursor: "pointer", pointerEvents: "auto" }}
      >
        {paths}
        {selection}
      </g>
    );
  }
}
