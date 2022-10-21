import React from "react";
import { connect } from "react-redux";
import { getAxisTitle } from "../reducers/plotData";
import { plotText } from "../reducers/plotStyles";
import styles from "./Plot.scss";

export default class AxisTitle extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => ({
        plotText: plotText(state),
      });
    };
  }

  render() {
    const ConnectedAxisTitle = connect(this.mapStateToProps)(
      AxisTitleComponent
    );
    return <ConnectedAxisTitle {...this.props} />;
  }
}

const AxisTitleComponent = ({
  axis,
  title,
  side = 1000,
  offset = 0,
  plotText,
}) => {
  let params = {};
  params.transform =
    axis == "x"
      ? `translate(${500 + offset},1085)`
      : `translate(-90,${500 - offset}),rotate(90)`;
  return (
    <g {...params}>
      <text style={plotText.axisTitle} transform={"scale(" + side / 1000 + ")"}>
        {title}
      </text>
    </g>
  );
};
