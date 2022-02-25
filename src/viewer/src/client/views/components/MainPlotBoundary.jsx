import { Axis, BOTTOM, LEFT, axisPropsFromTickScale } from "react-d3-axis";
import { plotPaths, plotText } from "../reducers/plotStyles";

import React from "react";
import { connect } from "react-redux";
import { format as d3Format } from "d3-format";
import { scaleLinear as d3scaleLinear } from "d3-scale";
import { getMainPlotData } from "../reducers/plotData";
import { getTransformFunctionParams } from "../reducers/plotParameters";

export default class MainPlotBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => ({
        data: getMainPlotData(state),
        params: getTransformFunctionParams(state),
        plotPaths: plotPaths(state),
        plotText: plotText(state),
      });
    };
  }

  render() {
    const PlotBoundary = connect(this.mapStateToProps)(PlotOutline);
    return <PlotBoundary {...this.props} />;
  }
}

const isPowerOfTen = (d) => {
  if (d > 1) {
    while (d > 9 && d % 10 == 0) {
      d /= 10;
    }
  } else if (d < 1) {
    return String(d).match(/^[0\.1]+$/);
  }
  return d == 1;
};

const PlotOutline = ({ data, params, plotPaths, plotText, fill }) => {
  let fontSize = plotText.axisTick.fontSize;
  let f = (d) => {
    if (d < 1 && d > 0.0001 && String(d).match(/^[0\.1]+$/)) {
      return d;
    }
    if (isPowerOfTen(d)) {
      return d3Format(".0s")(d);
    }
    return "";
  };
  let g = (d) => {
    if (d < 1 && d > 0.0001 && String(d).match(/^[0\.1]+$/)) {
      return d;
    }
    return d3Format(".0s")(d);
  };
  let xRange = data.meta.x.range;
  let xScale = data.meta.x.xScale.copy();
  let fx = f;
  if (data.meta.x.meta.scale != "scaleLog") {
    fx = (d) => d;
  } else if (xRange[0] * 100 > xRange[1]) {
    fx = g;
  }
  xScale.range([50, 950]);
  let xBreak, xBreakAxis;
  if (
    data.meta.x.meta.clamp &&
    data.meta.x.meta.clamp > data.meta.x.meta.limit[0]
  ) {
    let scale = d3scaleLinear()
      .range([50, 86])
      .domain([0, data.meta.x.meta.clamp]);
    xBreakAxis = (
      <g transform="translate(0,1010)">
        <Axis
          {...axisPropsFromTickScale(scale, 1)}
          style={{ orient: BOTTOM, tickFontSize: 0 }}
        />
        <text
          transform={"translate(" + scale(data.meta.x.meta.clamp * 0.5) + ",5)"}
          fontSize={fontSize}
          textAnchor="middle"
          dominantBaseline="hanging"
        >
          &lt; {data.meta.x.meta.clamp}
        </text>
      </g>
    );
    xBreak = (
      <line
        style={plotPaths.clampedDivider}
        x1={104}
        x2={104}
        y1={-300}
        y2={1050}
      />
    );
    xScale.range([122, 950]);
    xScale.domain([data.meta.x.meta.clamp, xScale.domain()[1]]);
  }
  let yRange = data.meta.y.range;
  let yScale = data.meta.y.xScale.copy();
  let fy = f;
  if (data.meta.y.meta.scale != "scaleLog") {
    fy = (d) => d;
  } else if (yRange[0] * 100 > yRange[1]) {
    fy = g;
  }
  yScale.range([950, 50]);
  let yBreak, yBreakAxis;
  if (
    data.meta.y.meta.clamp &&
    data.meta.y.meta.clamp > data.meta.y.meta.limit[0]
  ) {
    let scale = d3scaleLinear()
      .range([950, 914])
      .domain([0, data.meta.y.meta.clamp]);
    yBreakAxis = (
      <g transform="translate(-10)">
        <Axis
          {...axisPropsFromTickScale(scale, 1)}
          style={{ orient: LEFT, tickFontSize: 0 }}
        />
        <text
          transform={
            "translate(-5," + scale(data.meta.y.meta.clamp * 0.5) + ")"
          }
          fontSize={fontSize}
          textAnchor="end"
          dominantBaseline="middle"
        >
          &lt; {data.meta.y.meta.clamp}
        </text>
      </g>
    );
    yBreak = (
      <line
        style={plotPaths.clampedDivider}
        x1={-50}
        x2={1300}
        y1={896 + 50 * params.factor - params.intercept}
        y2={896 - 1300 * params.factor - params.intercept}
      />
    );
    yScale.range([878, 50]);
    yScale.domain([data.meta.y.meta.clamp, yScale.domain()[1]]);
  }

  return (
    <g>
      <rect
        style={{
          ...plotPaths.boundary,
          ...(fill && { fill }),
        }}
        x={0}
        y={0}
        width={1000}
        height={1000}
      />
      {xBreak}
      {xBreakAxis}
      {yBreak}
      {yBreakAxis}
      <Axis
        {...axisPropsFromTickScale(xScale, 10)}
        style={{ orient: BOTTOM, tickFontSize: 0 }}
      />
      <Axis
        {...axisPropsFromTickScale(yScale, 10)}
        style={{ orient: LEFT, tickFontSize: fontSize }}
        format={fy}
      />
      <g transform={"translate(1000)"}>
        <Axis
          {...axisPropsFromTickScale(yScale, 10)}
          style={{ orient: LEFT, tickFontSize: 0 }}
        />
      </g>
      <g transform={"translate(0,1000)"}>
        <Axis
          {...axisPropsFromTickScale(xScale, 10)}
          style={{ orient: BOTTOM, tickFontSize: fontSize }}
          format={fx}
        />
      </g>
    </g>
  );
};
