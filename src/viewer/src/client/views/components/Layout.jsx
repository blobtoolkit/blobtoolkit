import {
  getDatasetID,
  getHashString,
  getQueryString,
  getStatic,
  setQueryString,
} from "../reducers/location";
import {
  getDatasetIsActive,
  getReloading,
  getStaticThreshold,
  setReloading,
} from "../reducers/repository";

import BTKLogos from "./BTKLogos";
import CookieConsent from "react-cookie-consent";
import DOIBadge from "./DOIBadge";
import DatasetSpinner from "./DatasetSpinner";
import ExternalLink from "./ExternalLink";
import LayoutControls from "./LayoutControls";
import LayoutHeader from "./LayoutHeader";
import LayoutPlots from "./LayoutPlots";
import React from "react";
import Spinner from "./Spinner";
import { connect } from "react-redux";
import { gdprUrl } from "../reducers/gdpr";
import { getBinsForCat } from "../reducers/field";
import { getCatAxis } from "../reducers/plot";
import { getColorScheme } from "../reducers/color";
import { getSelectedDatasetMeta } from "../reducers/dataset";
import qs from "qs";
import { queryToStore } from "../querySync";
import styles from "./Layout.scss";

const branch = BRANCH || "";
const version = VERSION || "";
const git_version = GIT_VERSION || "";
const hash = COMMIT_HASH || "";

window.scrollTop = {};

class LayoutComponent extends React.Component {
  constructor(props) {
    super(props);
  }
  handleScroll(e, tab) {
    window.scrollTop[tab] = e.target.scrollTop;
  }

  render() {
    let message =
      "We use browser cookies to analyse traffic to this site. To use this site you must agree to our ";
    let url = "https://github.com/blobtoolkit/viewer";
    // {this.props.datasetId ? this.props.active ? <LayoutPlots/> : <Spinner/> : <LayoutPlots/> }
    let notice;
    if (gdprUrl) {
      notice = (
        <CookieConsent
          style={{ background: this.props.colors.darkColor }}
          contentStyle={{ margin: "15px" }}
          buttonText="Accept"
          declineButtonText="Decline"
          enableDeclineButton={true}
          cookieValue={true}
          declineCookieValue={false}
          buttonStyle={{
            backgroundColor: this.props.colors.highlightColor,
            color: this.props.colors.lightColor,
            border: `${this.props.colors.lightColor} solid 1px`,
            fontSize: "13px",
          }}
          declineButtonStyle={{
            backgroundColor: this.props.colors.lightColor,
            color: this.props.colors.highlightColor,
            border: `${this.props.colors.highlightColor} solid 1px`,
            fontSize: "13px",
          }}
          flipButtons={true}
        >
          This website uses cookies to help us monitor usage. To allow us to do
          this, please accept the terms of our{" "}
          <a
            style={{
              textDecoration: "underline",
              color: this.props.colors.highlightColor,
              fontWeight: "bold",
            }}
            href={gdprUrl}
            target="_blank"
          >
            Privacy Policy
          </a>
          .
        </CookieConsent>
      );
    }
    let plotsDiv;
    if (this.props.active != "loading") {
      plotsDiv = <LayoutPlots />;
    } else {
      plotsDiv = <LayoutPlots />;
      // plotsDiv = (
      //   <div className={styles.fill_parent}>
      //     <DatasetSpinner />
      //   </div>
      // );
    }

    return (
      <div className={styles.main}>
        <div className={styles.main_header}>
          <LayoutHeader />
        </div>
        <div className={styles.content}>
          <LayoutControls />
          <div className={styles.plot_area}>{plotsDiv}</div>
        </div>
        <div className={styles.main_footer}>
          <span style={{ float: "left" }}>
            BlobToolKit Viewer{" "}
            <a style={{ color: "white" }} href={url} target={git_version}>
              {version}
            </a>
          </span>
          <span style={{ float: "left", marginLeft: "1em" }}>
            <a
              style={{ color: "white" }}
              href="https://doi.org/10.1534/g3.119.400908"
              target="_blank"
            >
              Challis <i>et al.</i> 2020
            </a>
          </span>
          <BTKLogos />
          {notice}
        </div>
      </div>
    );
  }
}

class Layout extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      return {
        active: getDatasetIsActive(state),
        //queryString: getQueryString(state),
        reloading: getReloading(state),
        colors: getColorScheme(state),
        // bins: getBinsForCat(state),
        // cat: getCatAxis(state),
        // staticThreshold: getStaticThreshold(state),
        // meta: getSelectedDatasetMeta(state),
        // datasetId: getDatasetID(state),
        // isStatic: getStatic(state),
        //activeTab: getHashString(state)
      };
    };
    this.mapDispatchToProps = (dispatch) => {
      return {
        updateStore: (str, currentSearch, action) => {
          let values = qs.parse(str.replace("?", ""));
          dispatch(
            queryToStore({ values, searchReplace: true, currentSearch, action })
          );
        },
        updateQueryString: (qStr) => dispatch(setQueryString(qStr)),
      };
    };
  }

  render() {
    const ConnectedLayout = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(LayoutComponent);
    return <ConnectedLayout {...this.props} />;
  }
}

export default Layout;
