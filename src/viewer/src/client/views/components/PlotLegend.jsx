import {
  getCircleLimit,
  getLargeFonts,
  getPlotShape,
  getShowTotal,
} from "../reducers/plotParameters";
import { getDatasetID, getView } from "../reducers/location";

import React from "react";
import { connect } from "react-redux";
import { format as d3format } from "d3-format";
import { getDatasetMeta } from "../reducers/repository";
import { getSummary } from "../reducers/summary";
import { getWindowBinsForCat } from "../reducers/plotData";
import { plotText } from "../reducers/plotStyles";
import styles from "./Plot.scss";

export default class PlotLegend extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        let summary = getSummary(state);
        let id = getDatasetID(state);
        let meta = getDatasetMeta(state, id);
        let shape = getPlotShape(state);
        let view = getView(state);
        let circleLimit = getCircleLimit(state);
        let showTotal = getShowTotal(state);
        let largeFonts = getLargeFonts(state);
        return {
          meta,
          ...summary,
          shape,
          view,
          circleLimit,
          showTotal,
          largeFonts,
          bins: getWindowBinsForCat(state),
          plotText: plotText(state),
        };
      };
    };
  }

  render() {
    const ConnectedLegend = connect(this.mapStateToProps)(Legend);
    return <ConnectedLegend {...this.props} />;
  }
}

const Legend = ({
  values,
  zAxis,
  bins,
  palette,
  other,
  reducer,
  meta,
  shape,
  view,
  circleLimit,
  showTotal,
  largeFonts,
  plotText,
}) => {
  let items = [];
  let legendKey;
  let ds;
  let format = (x) => {
    if (x < 1000) {
      return x;
    }
    return d3format(".2s")(x);
  };
  let commaFormat = d3format(",");
  if (bins) {
    let offset = 20;
    let w = 19;
    let h = 19;
    let gap = 5;
    let title_x = w + gap;
    let title_y = offset - gap;
    if (largeFonts) {
      offset = 0;
      w = 21;
      h = 21;
      title_x = w + gap + 73;
      title_y = offset - 10;
    }

    ds = (
      <g transform={"translate(" + 0 + "," + 0 + ")"}>
        <text style={plotText.legendTitle}>{meta.name}</text>
      </g>
    );
    let headers = ["count"];
    if (reducer != "count") {
      headers.push(reducer + " " + zAxis);
    }
    if (zAxis == "length") {
      headers.push("n50");
    }
    legendKey = (
      <g transform={"translate(" + title_x + "," + title_y + ")"}>
        <text
          transform={"translate(" + (largeFonts ? -30 : 260) + ")"}
          style={Object.assign({}, plotText.legend, {
            textAnchor: largeFonts ? "start" : "end",
            fontWeight: "normal",
          })}
        >
          [{headers.join("; ")}]
        </text>
      </g>
    );
    let title = "total";
    let color = "#999";
    let numbers = [];
    let count = values.counts.all > 0;
    // numbers.push(commaFormat(values.counts.all))
    numbers.push(format(values.counts.all));
    if (reducer != "count") {
      numbers.push(format(values.reduced.all));
    }
    if (zAxis == "length") {
      numbers.push(format(values.n50.all));
    }
    if (count && showTotal) {
      items.push(
        <g
          key="all"
          transform={"translate(" + (largeFonts ? 101 : 0) + "," + offset + ")"}
        >
          {(largeFonts && (
            <g>
              <rect
                x={-155}
                y={-gap}
                width={150}
                height={h + gap}
                style={{ fill: "white", stroke: "none" }}
              />
              <text
                style={Object.assign({}, plotText.legend, {
                  textAnchor: "end",
                })}
                transform={"translate(" + -gap * 2 + "," + (h - gap) + ")"}
              >
                {title}
              </text>
              <text
                style={plotText.legend}
                transform={"translate(" + (w + gap * 2) + "," + (h - gap) + ")"}
              >
                [{numbers.join("; ")}]
              </text>
            </g>
          )) || (
            <g>
              <text
                style={plotText.legend}
                transform={"translate(" + (w + gap) + "," + (h - gap) + ")"}
              >
                {title}
              </text>
              <text
                style={Object.assign({}, plotText.legend, {
                  textAnchor: "end",
                })}
                transform={
                  "translate(" + (w + gap + 260) + "," + (h - gap) + ")"
                }
              >
                [{numbers.join("; ")}]
              </text>
            </g>
          )}
          <rect
            x={0}
            y={0}
            width={w}
            height={h}
            style={{ fill: color, stroke: "black" }}
          />
        </g>
      );
      offset += h + gap;
    }
    bins.forEach((bin, i) => {
      let title = bin.id;
      if (
        title == "no-hit" &&
        values.counts.all > circleLimit &&
        view == "blob" &&
        shape == "circle"
      ) {
        title += " (not shown)";
      }
      let color = palette.colors[i];
      let numbers = [];
      let count = values.counts.binned[i] > 0;
      // numbers.push(commaFormat(values.counts.binned[i]))
      numbers.push(format(values.counts.binned[i]));
      if (reducer != "count") {
        numbers.push(format(values.reduced.binned[i]));
      }
      if (zAxis == "length") {
        numbers.push(format(values.n50.binned[i]));
      }
      if (count || shape == "lines") {
        items.push(
          <g
            key={i}
            transform={
              "translate(" + (largeFonts ? 101 : 0) + "," + offset + ")"
            }
          >
            {(largeFonts && (
              <g>
                <rect
                  x={-155}
                  y={-gap}
                  width={325}
                  height={h + gap}
                  style={{ fill: "white", stroke: "none" }}
                />
                <text
                  style={Object.assign({}, plotText.legend, {
                    textAnchor: "end",
                  })}
                  transform={"translate(" + -gap * 2 + "," + (h - gap) + ")"}
                >
                  {title}
                </text>
                {count && (
                  <text
                    style={plotText.legend}
                    transform={
                      "translate(" + (w + gap * 2) + "," + (h - gap) + ")"
                    }
                  >
                    [{numbers.join("; ")}]
                  </text>
                )}
              </g>
            )) || (
              <g>
                <text
                  style={plotText.legend}
                  transform={"translate(" + (w + gap) + "," + (h - gap) + ")"}
                >
                  {title}
                </text>
                {count && (
                  <text
                    style={Object.assign({}, plotText.legend, {
                      textAnchor: "end",
                    })}
                    transform={
                      "translate(" + (w + gap + 260) + "," + (h - gap) + ")"
                    }
                  >
                    [{numbers.join("; ")}]
                  </text>
                )}
              </g>
            )}
            <rect
              x={0}
              y={0}
              width={w}
              height={h}
              style={{ fill: color, stroke: "black" }}
            />
          </g>
        );
        offset += h + gap;
      }
    });
  }
  return (
    <g>
      {largeFonts || ds}
      {legendKey}
      {items}
    </g>
  );
};
