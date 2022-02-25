import React, { memo } from "react";
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
  getDetailsForX,
  getDetailsForY,
  getScatterPlotData,
} from "../reducers/plotData";

import MainPlotBoundary from "./MainPlotBoundary";
import PlotParameters from "./PlotParameters";
import Pointable from "react-pointable";
import { connect } from "react-redux";
import { polygonContains as d3polygonContains } from "d3-polygon";
import { polygonHull as d3polygonHull } from "d3-polygon";
import { scaleLinear as d3scaleLinear } from "d3-scale";
import { getColorScheme } from "../reducers/color";
import { getPlotShape } from "../reducers/plotParameters";
import styles from "./Plot.scss";

export default class LassoSelect extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        return {
          data: getScatterPlotData(state),
          polygon: getSelectPolygon(state),
          selectedRecords: getSelectedRecords(state),
          selectionDisplay: getSelectionDisplay(state),
          xMeta: getDetailsForX(state),
          yMeta: getDetailsForY(state),
          selectSource: getSelectSource(state),
          colors: getColorScheme(state),
          shape: getPlotShape(state),
        };
      };
    };
    this.mapDispatchToProps = (dispatch) => {
      return {
        addRecords: (arr) => {
          dispatch(setSelectSource("circle"));
          dispatch(addRecords(arr));
        },
        replaceRecords: (arr) => {
          dispatch(setSelectSource("circle"));
          dispatch(replaceRecords(arr));
        },
        selectNone: (arr) => {
          dispatch(setSelectSource("circle"));
          dispatch(selectNone());
        },
        changeSelectPolygon: (arr) => {
          dispatch(setSelectSource("circle"));
          dispatch(changeSelectPolygon(arr));
        },
      };
    };
  }

  render() {
    const ConnectedLassoSelect = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(memo(Lasso));
    return <ConnectedLassoSelect {...this.props} />;
  }
}

const xScale = d3scaleLinear().range([-55.556, 955.556]);
const yScale = d3scaleLinear().range([-55.556, 955.556]);

const relativeCoords = (event) => {
  let bounds = event.target.getBoundingClientRect();
  //let width = bounds.right - bounds.left;
  xScale.domain([bounds.left, bounds.right]);
  yScale.domain([bounds.top, bounds.bottom]);
  let x = xScale(event.clientX);
  let y = yScale(event.clientY);
  //let height = bounds.bottom - bounds.top;
  //let x = ((event.clientX - bounds.left) - width * 0.05) / (width * 0.9) * 1000;
  //let y = ((event.clientY - bounds.top) - height * 0.05) / (height * 0.9) * 1000;
  return { x, y };
};

class Lasso extends React.Component {
  constructor(props) {
    super(props);

    let { polygon, internodes, complete } = this.make_polygon();
    this.state = {
      active: undefined,
      mouseDown: false,
      points: polygon,
      internodes,
      last: {},
      current: [],
      complete,
      moveNode: -1,
    };
  }

  make_polygon() {
    let polygon = this.props.polygon;
    let internodes = [];
    let complete = false;
    if (
      this.props.selectSource != "circle" &&
      this.props.selectSource != "list"
    ) {
      let len = this.props.selectedRecords.length;
      if (len > 2) {
        let points = [];
        for (let i = 0; i < len; i++) {
          let point = this.props.data.data.find((obj) => {
            return obj.id === this.props.selectedRecords[i];
          });
          if (point) points.push([point.x, point.y]);
        }
        polygon = d3polygonHull(points);
      } else {
        polygon = [];
      }
    }
    if (polygon && polygon.length > 2) {
      if (
        this.props.selectSource == "circle" ||
        this.props.selectSource == "list"
      ) {
        let xConvert = this.props.xMeta.xScale.copy().range([0, 900]);
        let yConvert = this.props.yMeta.xScale.copy().range([900, 0]);
        polygon = polygon.map((arr) => [xConvert(arr[0]), yConvert(arr[1])]);
      }
      if (this.props.selectSource == "list") {
        this.updateSelection(polygon, this.props.data.data);
      }
      complete = true;
      polygon.forEach((x, i) => {
        let point;
        if (i > 0) {
          point = [
            (x[0] + polygon[i - 1][0]) / 2,
            (x[1] + polygon[i - 1][1]) / 2,
          ];
        } else if (complete) {
          point = [
            (x[0] + polygon[polygon.length - 1][0]) / 2,
            (x[1] + polygon[polygon.length - 1][1]) / 2,
          ];
        }
        if (point) {
          internodes.push(point);
        }
      });
    }
    return { polygon, internodes, complete };
  }

  updateSelection(points, data) {
    let records = [];
    if (points.length > 2) {
      let len = data.length;
      for (let i = 0; i < len; i++) {
        let point = data[i];
        if (d3polygonContains(points, [point.x, point.y])) {
          records.push(point.id);
        }
      }
    }
    let xConvert = this.props.xMeta.xScale.copy().range([0, 900]);
    let yConvert = this.props.yMeta.xScale.copy().range([900, 0]);
    points = points.map((arr) => [
      xConvert.invert(arr[0]),
      yConvert.invert(arr[1]),
    ]);
    this.props.changeSelectPolygon(points);
    this.props.replaceRecords(records);
  }

  setMouseDown(bool, coords, marquee) {
    //this.setState({mouseDown:bool})
    let active = this.state.active;
    let arr = coords ? [coords.x, coords.y] : [];
    let complete = false;
    if (bool) {
      if (!active) {
        active = true;
        this.props.selectNone();
      }
    }
    if (bool && !marquee) {
      this.setState({
        mouseDown: true,
        current: arr,
        active: true,
      });
    } else {
      if (!coords) {
        if (this.state.complete) {
          return;
        }
        this.setState({
          mouseDown: bool,
          points: [],
          current: [],
          last: {},
          active: false,
          complete: false,
          internodes: [],
        });
        return;
      }
      let points = this.state.points.slice(0);
      let last = coords;
      if (this.state.hover >= 0) {
        if (this.state.hover == 0) {
          active = false;
          complete = true;
          if (this.state.clear) {
            clearTimeout(this.state.clear);
            this.setState({ clear: false });
          }
          // points.push(arr)
          this.updateSelection(points, this.props.data.data);
        } else {
          complete = false;
          last = {
            x: points[this.state.hover][0],
            y: points[this.state.hover][1],
          };
          points = points.slice(0, this.state.hover + 1);
        }
      } else {
        if (this.state.complete) {
          points = [arr];
          last = coords;
          this.props.selectNone();
        } else {
          points.push(arr);
        }
      }
      // if (this.state.points.length > 1){
      //   points = this.state.points.concat([this.state.points[0]])
      //   last = {
      //     x:this.state.points[this.state.points.length-1][0],
      //     y:this.state.points[this.state.points.length-1][1]}
      // }
      let internodes = [];
      points.forEach((x, i) => {
        let point;
        if (i > 0) {
          point = [
            (x[0] + points[i - 1][0]) / 2,
            (x[1] + points[i - 1][1]) / 2,
          ];
        } else if (complete) {
          point = [
            (x[0] + points[points.length - 1][0]) / 2,
            (x[1] + points[points.length - 1][1]) / 2,
          ];
        }
        if (point) {
          internodes.push(point);
        }
      });
      this.setState({
        mouseDown: bool,
        points,
        current: [],
        last,
        active,
        complete,
        internodes,
      });
    }
  }

  clearIncomplete() {
    if (this.state.clear) {
      clearTimeout(this.state.clear);
    }
    if (this.state.complete) {
      return;
    }

    let clear = setTimeout(() => {
      this.setState({
        mouseDown: false,
        points: [],
        current: [],
        last: {},
        active: false,
        complete: false,
        internodes: [],
      });
    }, 3000);
    this.setState({ clear });
  }

  componentDidUpdate(nextProps) {
    if (
      nextProps.selectSource != "circle" &&
      nextProps.selectedRecords != this.props.selectedRecords
    ) {
      let { polygon, internodes, complete } = this.make_polygon();
      this.setState({ points: polygon, internodes, complete });
    }
  }

  componentWillUnmount() {
    // fix Warning: Can't perform a React state update on an unmounted component
    this.setState = (state, callback) => {
      return;
    };
  }

  drawNode(coords, key, fill = -1) {
    let nodeWidth = 20;
    return (
      <rect
        fill={
          fill ? this.props.colors.highlightColor : this.props.colors.lightColor
        }
        key={key}
        x={coords[0] - nodeWidth / 2}
        y={coords[1] - nodeWidth / 2}
        height={nodeWidth}
        width={nodeWidth}
      />
    );
  }

  drawInterNode(coords, key, fill = -1) {
    let nodeWidth = 14;
    return (
      <rect
        fill={
          fill ? this.props.colors.highlightColor : this.props.colors.lightColor
        }
        key={key}
        x={coords[0] - nodeWidth / 2}
        y={coords[1] - nodeWidth / 2}
        fillOpacity={fill ? 1 : 0.5}
        height={nodeWidth}
        width={nodeWidth}
      />
    );
  }

  detectNodes(coords, thresh) {
    let hover = -1;
    let ihover = -1;
    let current = [coords.x, coords.y];
    this.state.points.forEach((point, i) => {
      let deltaX = coords.x - point[0];
      let deltaY = coords.y - point[1];
      if (Math.abs(deltaX) <= thresh && Math.abs(deltaY) <= thresh) {
        current = [point[0], point[1]];
        hover = i;
      }
    });
    this.state.internodes.forEach((point, i) => {
      let deltaX = coords.x - point[0];
      let deltaY = coords.y - point[1];
      if (Math.abs(deltaX) <= thresh && Math.abs(deltaY) <= thresh) {
        ihover = i;
      }
    });
    let state = {
      current,
      hover,
      ihover,
    };
    this.setState(state);
    return state;
  }

  render() {
    let points = this.state.points;
    let thresh = 15;
    let latest;
    let nodes, internodes;
    let path;

    if (
      this.props.selectionDisplay // &&
      // this.props.selectSource == this.props.shape
    ) {
      if (this.state.active) {
        let current = this.state.current;
        let line;
        if (current.length == 2) {
          if (points.length > 0) {
            line = (
              <path
                d={`M${this.state.last.x} ${this.state.last.y}L${current.join(
                  " "
                )}`}
              />
            );
          }
          let square = this.drawNode(current, 0, false);

          latest = (
            <g
              fill="none"
              stroke={this.props.colors.highlightColor}
              strokeWidth={5}
              opacity={0.5}
            >
              {line}
              {square}
            </g>
          );
        }
      }
      if (points && points.length > 0) {
        nodes = [];
        points.forEach((x, i) => {
          nodes.push(this.drawNode(x, i, this.state.hover == i));
        });
        internodes = [];
        this.state.internodes.forEach((x, i) => {
          internodes.push(this.drawInterNode(x, i, this.state.ihover == i));
        });
        if (points.length > 1) {
          let complete = this.state.complete ? "Z" : "";
          path = (
            <path
              stroke={this.props.colors.highlightColor}
              strokeWidth={5}
              d={`M${(points[0] || []).join(" ")} ${(points.slice(1) || [])
                .map((x) => "L" + (x || []).join(" "))
                .join("")}${complete}`}
            />
          );
        }
      }
    }
    let dasharray = this.props.selectSource == "circle" ? undefined : "10";
    let fill = this.state.complete
      ? this.props.selectSource == "circle"
        ? this.props.colors.halfHighlightColor
        : "url(#pattern-checkers)"
      : "none";
    return (
      <g>
        <pattern
          id="pattern-checkers"
          x="0"
          y="0"
          width="20"
          height="20"
          patternUnits="userSpaceOnUse"
        >
          <rect
            x="0"
            width="10"
            height="10"
            y="0"
            fill={this.props.colors.halfHighlightColor}
          />
          <rect
            x="10"
            width="10"
            height="10"
            y="10"
            fill={this.props.colors.halfHighlightColor}
          />
        </pattern>
        {this.props.selectionDisplay && (
          <g transform="translate(50,50)">
            {latest}
            {points && points.length > 1 && (
              <g strokeDasharray={dasharray} fill={fill}>
                {path}
              </g>
            )}
            {points && points.length > 0 && (
              <g
                fill={this.props.colors.lightColor}
                stroke={this.props.colors.highlightColor}
                strokeWidth={5}
              >
                {nodes}
                {internodes}
              </g>
            )}
          </g>
        )}
        <Pointable
          tagName="g"
          onPointerMove={(e) => {
            e.preventDefault();
            let coords = relativeCoords(e);
            let current = [coords.x, coords.y];
            if (this.state.moveNode >= 0) {
              points[this.state.moveNode] = current;
              let internodes = [];
              points.forEach((x, i) => {
                let point;
                if (i > 0) {
                  point = [
                    (x[0] + points[i - 1][0]) / 2,
                    (x[1] + points[i - 1][1]) / 2,
                  ];
                } else if (this.state.complete) {
                  point = [
                    (x[0] + points[points.length - 1][0]) / 2,
                    (x[1] + points[points.length - 1][1]) / 2,
                  ];
                }
                if (point) {
                  internodes.push(point);
                }
              });
              this.setState({
                points,
                internodes,
              });
            } else {
              // if (this.state.mouseDown) {
              //   if (points.length == 0) {
              //     this.setMouseDown(true, coords, true);
              //   }
              //   if (points.length >= 1) {
              //     let newPoints = [
              //       points[0],
              //       [current[0], points[0][1]],
              //       current,
              //       [points[0][0], current[1]],
              //       points[0],
              //     ];
              //     this.setState({
              //       points: newPoints,
              //     });
              //   }
              // } else {
              this.detectNodes(coords, thresh);
              this.clearIncomplete();
              // }
            }
          }}
          onPointerLeave={(e) => {
            e.preventDefault();
            this.setMouseDown(false);
          }}
          onPointerDown={(e) => {
            e.preventDefault();
            let coords = relativeCoords(e);
            let state = this.detectNodes(coords, thresh);
            if (this.state.hover >= 0 || state.hover >= 0) {
              let hover = Math.max(this.state.hover, state.hover);
              if (this.state.complete) {
                this.setState({
                  hover,
                  moveNode: hover,
                });
                return;
              } else {
                this.setState({
                  hover,
                });
              }
            } else if (
              this.state.complete &&
              (this.state.ihover >= 0 || state.ihover >= 0)
            ) {
              let ihover = Math.max(this.state.ihover, state.ihover);
              let hover = this.state.ihover;
              let moveNode = hover;
              ihover = -1;
              points.splice(
                this.state.ihover,
                0,
                this.state.internodes[this.state.ihover]
              );
              this.setState({
                points,
                hover,
                moveNode,
                ihover,
              });
              return;
            }
            this.setMouseDown(true, coords);
            this.clearIncomplete();
          }}
          onPointerUp={(e) => {
            e.preventDefault();
            if (this.state.moveNode >= 0) {
              this.updateSelection(this.state.points, this.props.data.data);
              this.setState({ moveNode: -1 });
            } else {
              let coords = relativeCoords(e);
              if (this.state.hover >= 0) {
                coords = {
                  x: this.state.points[this.state.hover][0],
                  y: this.state.points[this.state.hover][1],
                };
              }
              this.setMouseDown(false, coords);
              this.clearIncomplete();
            }
          }}
        >
          <MainPlotBoundary />
        </Pointable>
      </g>
    );
  }
}
