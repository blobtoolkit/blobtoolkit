import { chooseView, getStatic, getView } from "../reducers/location";

import MenuDisplaySimple from "./MenuDisplaySimple";
import React from "react";
import TextIcon from "./TextIcon";
import { connect } from "react-redux";
import { getMainPlotData } from "../reducers/plotData";
import styles from "./Layout.scss";

class SwitchView extends React.Component {
  constructor(props) {
    super(props);
  }

  render() {
    let { view, isStatic, data, onSelectView } = this.props;
    if (!data) return null;
    let showBlob;
    if (data.axes && data.axes.y.values.length > 0) {
      showBlob = true;
    }
    return (
      <MenuDisplaySimple invert={true}>
        {showBlob && (
          <TextIcon
            invert={true}
            title="blob"
            active={view == "blob"}
            onIconClick={() => onSelectView("blob")}
          />
        )}
        <TextIcon
          invert={true}
          title="busco"
          active={view == "busco"}
          onIconClick={() => onSelectView("busco")}
        />
        <TextIcon
          invert={true}
          title="cumulative"
          active={view == "cumulative"}
          onIconClick={() => onSelectView("cumulative")}
        />
        <TextIcon
          invert={true}
          title="detail"
          active={view == "detail"}
          onIconClick={() => onSelectView("detail")}
        />
        {isStatic || (
          <TextIcon
            invert={true}
            title="report"
            active={view == "report"}
            onIconClick={() => onSelectView("report")}
          />
        )}
        <TextIcon
          invert={true}
          title="snail"
          active={view == "snail"}
          onIconClick={() => onSelectView("snail")}
        />
        {isStatic || (
          <TextIcon
            invert={true}
            title="table"
            active={view == "table"}
            onIconClick={() => onSelectView("table")}
          />
        )}
      </MenuDisplaySimple>
    );
  }
}

class MenuSwitchView extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      return {
        view: getView(state),
        isStatic: getStatic(state),
        data: getMainPlotData(state),
      };
    };
    this.mapDispatchToProps = (dispatch) => {
      return {
        onSelectView: (view) => dispatch(chooseView(view)),
      };
    };
  }

  render() {
    const ConnectedSwitchView = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(SwitchView);
    return <ConnectedSwitchView {...this.props} />;
  }
}

export default MenuSwitchView;
