import { Axis, BOTTOM, LEFT, axisPropsFromTickScale } from "react-d3-axis";
import { plotPaths, plotText } from "../reducers/plotStyles";

import React from "react";
import { connect } from "react-redux";
import { format as d3Format } from "d3-format";
import { scaleLinear as d3scaleLinear } from "d3-scale";
import { getMainPlotData } from "../reducers/plotData";
import { getSelectionDisplay } from "../reducers/select";
import { getTransformFunctionParams } from "../reducers/plotParameters";

export default class PlotGridBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => ({
        plotText: plotText(state),
        plotPaths: plotPaths(state),
        showSelection: getSelectionDisplay(state),
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

const PlotOutline = ({
  i,
  col,
  row,
  x,
  y,
  width,
  height,
  label,
  colWidth,
  xDomain,
  yDomain,
  scales,
  nRows,
  nCols,
  len,
  plotText,
  plotPaths,
  fill = "white",
  selected,
  handleClick,
  highlightColor,
  showSelection,
  meta,
}) => {
  let fontSize = plotText.axisTick.fontSize;
  let f = (d) => {
    if (d < 1 && d > 0.0001) {
      // && String(d).match(/^[0\.1]+$/)) {
      return String(d).replace(/\.0+$/, "");
    }
    if (isPowerOfTen(d)) {
      return d3Format(".2s")(d);
    }
    return "";
  };
  let g = (d) => {
    if (d < 1 && d > 0.0001) {
      return String(d).replace(/\.0+$/, "");
    }
    return d3Format(".2s")(d).replace(/\.0+$/, "");
  };
  let xScale = scales.x.copy();
  let fx = g;
  let xRange = xScale.domain();
  if (meta.x.meta.scale == "scaleLog" && xRange[0] * 100 <= xRange[1]) {
    fx = f;
  }
  xScale.range([0, width]);
  xScale.domain([xScale.domain()[0], colWidth]);
  let xBreak, xBreakAxis;
  // if (
  //   data.meta.x.meta.clamp &&
  //   data.meta.x.meta.clamp > data.meta.x.meta.limit[0]
  // ) {
  //   let scale = d3scaleLinear()
  //     .range([50, 86])
  //     .domain([0, data.meta.x.meta.clamp]);
  //   xBreakAxis = (
  //     <g transform="translate(0,1010)">
  //       <Axis
  //         {...axisPropsFromTickScale(scale, 1)}
  //         style={{ orient: BOTTOM, tickFontSize: 0 }}
  //       />
  //       <text
  //         transform={"translate(" + scale(data.meta.x.meta.clamp * 0.5) + ",5)"}
  //         fontSize={fontSize}
  //         textAnchor="middle"
  //         dominantBaseline="hanging"
  //       >
  //         &lt; {data.meta.x.meta.clamp}
  //       </text>
  //     </g>
  //   );
  //   xBreak = (
  //     <line
  //       style={plotPaths.clampedDivider}
  //       x1={104}
  //       x2={104}
  //       y1={-300}
  //       y2={1050}
  //     />
  //   );
  //   xScale.range([122, 950]);
  //   xScale.domain([data.meta.x.meta.clamp, xScale.domain()[1]]);
  // }
  let yScale = scales.y.copy();
  let fy = g;

  let yRange = yScale.domain();
  if (meta.y.meta.scale == "scaleLog" && yRange[0] * 100 <= yRange[1]) {
    fy = f;
  }
  yScale.range([height, 0]);
  let yBreak, yBreakAxis;
  // if (
  //   data.meta.y.meta.clamp &&
  //   data.meta.y.meta.clamp > data.meta.y.meta.limit[0]
  // ) {
  //   let scale = d3scaleLinear()
  //     .range([950, 914])
  //     .domain([0, data.meta.y.meta.clamp]);
  //   yBreakAxis = (
  //     <g transform="translate(-10)">
  //       <Axis
  //         {...axisPropsFromTickScale(scale, 1)}
  //         style={{ orient: LEFT, tickFontSize: 0 }}
  //       />
  //       <text
  //         transform={
  //           "translate(-5," + scale(data.meta.y.meta.clamp * 0.5) + ")"
  //         }
  //         fontSize={fontSize}
  //         textAnchor="end"
  //         dominantBaseline="middle"
  //       >
  //         &lt; {data.meta.y.meta.clamp}
  //       </text>
  //     </g>
  //   );
  //   // yBreak = (
  //   //   <line
  //   //     style={plotPaths.clampedDivider}
  //   //     x1={-50}
  //   //     x2={1300}
  //   //     y1={896 + 50 * params.factor - params.intercept}
  //   //     y2={896 - 1300 * params.factor - params.intercept}
  //   //   />
  //   // );
  //   yScale.range([878, 50]);
  //   yScale.domain([data.meta.y.meta.clamp, yScale.domain()[1]]);
  // }

  let xTicks = nCols == 1 ? 10 : nCols == 2 ? 5 : 3;
  let yTicks = nRows == 1 ? 10 : nRows == 2 ? 5 : 3;

  return (
    <g transform={`translate(${x},${y - height})`}>
      <rect
        // style={{
        //   // ...plotPaths.boundary,
        //   ...(fill && { fill }),
        // }}
        stroke={showSelection && selected ? highlightColor : "none"}
        strokeWidth={showSelection && selected ? "2px" : 0}
        fill={showSelection && selected ? highlightColor : "white"}
        fillOpacity={showSelection && selected ? 0.25 : 1}
        x={0}
        y={0}
        width={width}
        height={height}
        onPointerDown={handleClick}
      />
      <g transform={`translate(0,${10})`} pointerEvents={"none"}>
        <text
          x={width / 2}
          textAnchor={"middle"}
          dominantBaseline={"hanging"}
          alignmentBaseline={"hanging"}
          fontSize={fontSize}
        >
          {label}
        </text>
      </g>
      {/* {xBreak}
      {xBreakAxis}
      {yBreak}
      {yBreakAxis} */}
      {/* <Axis
        {...axisPropsFromTickScale(xScale, xTicks)}
        style={{ orient: BOTTOM, tickFontSize: 0 }}
      /> */}
      {(col == 0 && (
        <Axis
          {...axisPropsFromTickScale(yScale, yTicks)}
          style={{ orient: LEFT, tickFontSize: fontSize }}
          format={fy}
        />
      )) || (
        <Axis
          {...axisPropsFromTickScale(yScale, yTicks)}
          style={{ orient: LEFT, tickFontSize: 0 }}
        />
      )}
      {/* <g transform={`translate(${width})`}>
        <Axis
          {...axisPropsFromTickScale(yScale, yTicks)}
          style={{ orient: LEFT, tickFontSize: 0 }}
        />
      </g> */}
      {((row == nRows - 1 || i == len - 1) && (
        <g transform={`translate(0,${height})`}>
          <Axis
            {...axisPropsFromTickScale(xScale, xTicks)}
            style={{ orient: BOTTOM, tickFontSize: fontSize }}
            format={fx}
          />
        </g>
      )) || (
        <g transform={`translate(0,${height})`}>
          <Axis
            {...axisPropsFromTickScale(xScale, xTicks)}
            style={{ orient: BOTTOM, tickFontSize: 0 }}
          />
        </g>
      )}
    </g>
  );
};
