import { getDatasetID, getHashString } from "../reducers/location";

import MenuDatasetMain from "./MenuDatasetMain";
import MenuDisplayMain from "./MenuDisplayMain";
import MenuFilterMain from "./MenuFilterMain";
import MenuHelpMain from "./MenuHelpMain";
import MenuListsMain from "./MenuListsMain";
import MenuSummaryMain from "./MenuSummaryMain";
import React from "react";
import { connect } from "react-redux";
import { datasetTable } from "../reducers/api";
import { getTopLevelFields } from "../reducers/field";
import styles from "./Layout.scss";

class ControlsLayoutComponent extends React.Component {
  constructor(props) {
    super(props);
  }
  handleScroll() {
    let menuDiv = this.refs.menuDiv;
    let activeTab = this.props.activeTab;
    window.scrollTop[activeTab] = menuDiv.scrollTop;
  }

  componentDidMount() {
    let menuDiv = this.refs.menuDiv;
    if (menuDiv) {
      let activeTab = this.props.activeTab;
      menuDiv.scrollTop = window.scrollTop[activeTab] || 0;
      menuDiv.addEventListener("scroll", () => this.handleScroll());
    }
  }

  componentDidUpdate() {
    let menuDiv = this.refs.menuDiv;
    if (menuDiv) {
      let activeTab = this.props.activeTab;
      menuDiv.removeEventListener("scroll", () => this.handleScroll());
      menuDiv.scrollTop = window.scrollTop[activeTab] || 0;
      menuDiv.addEventListener("scroll", () => this.handleScroll());
    }
  }

  render() {
    let labels = [
      "Datasets",
      "Filters",
      "Lists",
      "Settings",
      "Summary",
      "Help",
    ];
    let tabs = [];
    let activeTab = this.props.activeTab;
    if (!activeTab && !this.props.datasetId) {
      activeTab = "Datasets";
    }
    labels.forEach((tab) => {
      tabs.push({ id: tab, active: activeTab == tab });
    });
    let menu;
    if (!datasetTable) {
      if (activeTab == "Datasets") menu = <MenuDatasetMain />;
    }
    if (activeTab == "Filters") menu = <MenuFilterMain />;
    if (activeTab == "Lists") menu = <MenuListsMain />;
    if (activeTab == "Settings") menu = <MenuDisplayMain />;
    if (activeTab == "Summary") menu = <MenuSummaryMain />;
    if (activeTab == "Help") menu = <MenuHelpMain />;
    if (menu) {
      return (
        <div className={styles.menu} ref="menuDiv">
          {menu}
        </div>
      );
    } else {
      return null;
    }
  }
}

class LayoutControls extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      return {
        topLevelFields: getTopLevelFields(state),
        activeTab: getHashString(state),
        datasetId: getDatasetID(state),
      };
    };
  }

  render() {
    const ConnectedLayout = connect(this.mapStateToProps)(
      ControlsLayoutComponent
    );
    return <ConnectedLayout {...this.props} />;
  }
}

export default LayoutControls;
