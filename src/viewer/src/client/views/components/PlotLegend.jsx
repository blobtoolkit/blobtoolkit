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

export const charWidth = (char, options = { factor: 0.7 }) => {
  const { factor } = options;
  const widths = {
    dot: 2,
    number: 5,
    a: 5,
    g: 7,
    i: 4,
    j: 4.5,
    m: 8,
    M: 10,
  };
  const chars = {};
  [".", ",", ";", ":", "|", "!", "\\", "/", " "].forEach((char) => {
    chars[char] = widths.dot;
  });
  [...Array(10).keys()].forEach((char) => {
    chars[char] = widths.number;
  });
  ["g"].forEach((char) => {
    chars[char] = widths.g;
  });
  ["i", "l", 1].forEach((char) => {
    chars[char] = widths.i;
  });
  ["j", "t", "-"].forEach((char) => {
    chars[char] = widths.j;
  });
  ["m", "w"].forEach((char) => {
    chars[char] = widths.m;
  });
  ["M", "W"].forEach((char) => {
    chars[char] = widths.M;
  });
  let width = widths.a;
  if (chars[char]) {
    width = chars[char];
  } else if (chars[char.toLowerCase()]) {
    width = widths.m;
  }
  return (width / widths.a) * factor;
};

const stringLength = (str, options) => {
  let length = `${isNaN(str) ? (str ? str : "") : str}`
    .split("")
    .reduce((a, b) => a + charWidth(b, options), 0);
  return length;
};

const Legend = ({
  values,
  zAxis,
  bins,
  visibleCats,
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
  compactLegend,
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
    let offsetY = gap;
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
    let legendRows = 0;
    // numbers.push(commaFormat(values.counts.all))
    if (shape == "grid") {
      offset = 0;
      let visibleBins = bins.filter((b) =>
        visibleCats ? visibleCats.has(b.id) : true
      );
      let fontSize = plotText.legend.fontSize.replace("px", "") * 1.2;
      let length = Math.max(
        ...visibleBins.map(
          (b) => stringLength(b.id, { factor: 0.6 }) * fontSize
        )
      );
      legendRows = Math.ceil(((length + w + gap) * bins.length) / 1200);
      if (legendRows > 3) {
        fontSize /= 1.2;
      }
      bins.forEach((bin, i) => {
        let title = bin.id;
        if (visibleCats.has(title)) {
          let color = palette.colors[i];
          if (legendRows > 2) {
            length = stringLength(title, { factor: 0.65 }) * fontSize;
          }
          if (offset + length > 1200) {
            offset = 0;
            offsetY += h + gap;
          }
          items.push(
            <g
              key={title}
              transform={"translate(" + offset + "," + offsetY + ")"}
            >
              <text
                style={Object.assign({}, plotText.legend, {
                  textAnchor: "start",
                  fontSize: `${fontSize}px`,
                  alignmentBaseline: "middle",
                  dominantBaseline: "middle",
                })}
                // fontSize={fontSize * 1.2}
                // alignmentBaseline={"middle"}
                // dominantBaseline={"middle"}
                transform={"translate(" + (w + gap) + "," + h / 2 + ")"}
              >
                {title}
              </text>

              <rect
                x={0}
                y={0}
                width={w}
                height={h}
                style={{ fill: color, stroke: "black" }}
              />
            </g>
          );
          offset += length + w + gap;
        }
      });
      return (
        <g transform={`translate(0,${offsetY == gap ? h / 2 : 0})`}>{items}</g>
      );
    }
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
