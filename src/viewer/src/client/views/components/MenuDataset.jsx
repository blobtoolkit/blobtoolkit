import { DatasetModal } from "./DatasetModal";
import { Link } from "react-router-dom";
import PropTypes from "prop-types";
import React from "react";
import Spinner from "./Spinner";
import { connect } from "react-redux";
import { createSelector } from "reselect";
import { getDatasetMeta } from "../reducers/repository";
import styles from "./Layout.scss";

const basename = BASENAME || "";

class MenuDataset extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => ({
      meta: getDatasetMeta(state, this.props.id),
    });
  }

  componentDidMount() {
    this.props.onDatasetMount(this.props.id);
  }

  render() {
    const ConnectedDataset = connect(this.mapStateToProps)(Dataset);
    return <ConnectedDataset {...this.props} />;
  }
}

class Dataset extends React.Component {
  constructor(props) {
    super(props);
    this.state = { show: false };
  }
  render() {
    if (!this.props.meta) {
      return <Spinner />;
    }
    let css = styles.menu_item;
    if (this.props.active) css += " " + styles.active;
    let records = this.props.meta.records || "";
    let record_type = this.props.meta.record_type || "";
    if (record_type && record_type.substr(-1) != "s") record_type += "s";
    let taxon_meta = this.props.meta.taxon || {};
    let taxon = taxon_meta.name || "";
    let accession = this.props.meta.accession || "";
    return (
      <div className={css}>
        <a
          data-tip
          data-for="load-dataset"
          className={styles.cover}
          onClick={() => this.props.onDatasetClick(this.props.id)}
        >
          <h3>{this.props.meta.name}</h3>
          <span className={styles.menu_subtitle}>
            {accession}
            <br />
            <em>{taxon}</em>
          </span>
        </a>
        <span className={styles.menu_count}>
          {records.toLocaleString() + " " + record_type}
        </span>
        <div
          data-tip
          data-for="view-metadata"
          className={styles.right}
          onClick={() => this.setState({ show: true })}
        >
          details
          <DatasetModal
            meta={this.props.meta}
            selected={this.state.show}
            dismiss={() => this.setState({ show: false })}
          />
        </div>
      </div>
    );
  }
}

export default MenuDataset;
