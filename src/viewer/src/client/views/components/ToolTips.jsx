import React from "react";
import ReactTooltip from "react-tooltip";
import { connect } from "react-redux";
import { editFilter } from "../reducers/filter";
import { getDetailsForFilterId } from "../reducers/preview";
import isTouchDevice from "is-touch-device";

const sets = (set) => {
  let tips = {};

  tips.xAxis = {
    name: "xAxis",
    detail: "assign variable to x-axis",
    filterMenu: 1,
  };
  tips.yAxis = {
    name: "yAxis",
    detail: "assign variable to y-axis",
    filterMenu: 1,
  };
  tips.zAxis = {
    name: "zAxis",
    detail: "use variable to scale points/bins",
    filterMenu: 1,
  };
  tips.category = { name: "category", filterMenu: 1 };
  tips.clone = { name: "sum active fields", filterMenu: 1 };
  tips.scale = { name: "toggle scale", filterMenu: 1 };
  tips.zero = { name: "curves start at 0", filterMenu: 1 };
  tips.xLetter = { name: "stack curves by x position", filterMenu: 1 };
  tips.yLetter = { name: "stack curves by y position", filterMenu: 1 };
  tips.sumActive = { name: "sum active fields", filterMenu: 1 };
  tips.invert = { name: "invert", filterMenu: 1 };
  tips.show = { name: "showSelection", filterMenu: 1 };
  tips.showHide = { name: "hideSelection", filterMenu: 1 };
  tips.selectAll = { name: "selectAll", filterMenu: 1 };
  tips.selectNone = { name: "selectNone", filterMenu: 1 };
  tips.invertSelection = { name: "invertSelection", filterMenu: 1 };
  tips["field-header"] = { name: "activate", filterMenu: 1 };
  //tips['header-type'] = {name:'field type',filterMenu:1}
  tips["active-field-header"] = { name: "deactivate", filterMenu: 1 };
  tips["draggable-arrow"] = { name: "drag to set range", filterMenu: 1 };
  tips["range-input"] = { name: "set range", filterMenu: 1 };
  tips["clamp-input"] = { name: "set clamp", filterMenu: 1 };
  tips["create-list-input"] = { name: "set list name", filterMenu: 1 };
  tips["create-list-button"] = {
    name: "convert current filter to list",
    filterMenu: 1,
  };
  tips["filter-summary"] = {
    name: "summary of filtered dataset",
    filterMenu: 1,
  };
  tips["filter-button"] = { name: "apply filter", filterMenu: 1 };
  tips["category-toggle"] = { name: "toggle category", filterMenu: 1 };

  tips.squareShape = { name: "plot square bins", settingsMenu: 1 };
  tips.hexShape = { name: "plot hexagonal bins", settingsMenu: 1 };
  tips.circleShape = { name: "plot points as circles", settingsMenu: 1 };
  tips.lineShape = { name: "plot connected windows", settingsMenu: 1 };
  tips.kiteShape = { name: "plot kite shapes", settingsMenu: 1 };
  tips.max = { name: "maximum z-value", settingsMenu: 1 };
  tips.min = { name: "minimum z-value", settingsMenu: 1 };
  tips.sum = { name: "sum of z-values", settingsMenu: 1 };
  tips.count = { name: "count", settingsMenu: 1 };
  tips.mean = { name: "mean z-value", settingsMenu: 1 };
  tips.log = { name: "log scale", settingsMenu: 1 };
  tips.linear = { name: "linear scale", settingsMenu: 1 };
  tips.adjust = { name: "exclude Ns from coverage", settingsMenu: 1 };
  tips.ellipsis = { name: "list reference assemblies", settingsMenu: 1 };
  tips.fontSize = { name: "toggle larger fonts", settingsMenu: 1 };
  tips.sqrt = { name: "squre-root scale", settingsMenu: 1 };
  tips["size-slider"] = { name: "scale factor for bin size", filterMenu: 1 };
  tips["scale-slider"] = {
    name: "scale factor for size within bin",
    filterMenu: 1,
  };
  tips["select-palette"] = { name: "select palette", filterMenu: 1 };
  tips["edit-swatch"] = { name: "adjust colour", filterMenu: 1 };
  tips["load-dataset"] = { name: "load dataset", datasetMenu: 1 };
  tips["view-metadata"] = { name: "view metadata", datasetMenu: 1 };

  return tips;
};

const ActiveToolTips = ({ set, detail }) => {
  let tips = sets(set);
  let toolTips = Object.keys(tips).map((tip) => (
    <ReactTooltip key={tip} id={tip}>
      {tips[tip].name}
    </ReactTooltip>
  ));
  return <span>{isTouchDevice() ? "" : toolTips}</span>;
};

class ToolTips extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = () => {
      return (state, props) => {
        return {};
        // return getDetailsForFilterId(state, props.filterId)
      };
    };
    this.mapDispatchToProps = (dispatch) => {
      return {
        // onUpdateRange: (id,range) => {return dispatch(editFilter({id,range}))}
      };
    };
  }

  render() {
    const ConnectedToolTips = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(ActiveToolTips);
    return <ConnectedToolTips {...this.props} />;
  }
}

export default ToolTips;
