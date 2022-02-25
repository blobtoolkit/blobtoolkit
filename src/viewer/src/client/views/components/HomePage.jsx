import FindDatasets from "./FindDatasets";
import React from "react";
import Search from "./Search";
import { connect } from "react-redux";
import figure1 from "./img/figure1.jpg";
import figure2 from "./img/figure2.jpg";
import hexBlur from "./img/hex_blur.jpg";
import styles from "./HomePage.scss";

const dataset_table = DATASET_TABLE || false;

export default class HomePage extends React.Component {
  constructor(props) {
    super(props);
  }

  render() {
    return (
      <div className={styles.outer}>
        <div
          style={{
            backgroundImage: `url(${hexBlur})`,
            margin: "-1em -2em 0 -2em",
            padding: "1em 2em 1em 2em",
            color: "white",
          }}
        >
          <h1>BlobToolKit Viewer</h1>
          <p>Interactive visualisation of genomic datasets.</p>
        </div>
        {dataset_table && <FindDatasets />}
        {!dataset_table && (
          <span className={styles.no_table}>
            Search datasets to begin (search 'all' to show all available
            datasets).
            <hr />
            <p>
              Multiple views and data export options dynamically update as
              filter parameters and selections are modified (Figure 1).
            </p>
            <img src={figure1} alt="Figure 1" />
            <p>
              We are running the{" "}
              <a href="https://github.com/blobtoolkit/insdc-pipeline">
                BlobToolKit pipeline
              </a>{" "}
              (Figure 2) on all public (INSDC registered) eukaryote genome
              assemblies and making the results available at{" "}
              <a href="http://blobtoolkit.genomehubs.org/view">
                blobtoolkit.genomehubs.org/view
              </a>
              .
            </p>
            <img src={figure2} alt="Figure 2" />
            <p>
              To find out more about BlobToolKit, visit the project homepage at{" "}
              <a href="http://blobtoolkit.genomehubs.org">
                blobtoolkit.genomehubs.org
              </a>{" "}
              or browse the FAQs in the{" "}
              <span
                className={styles.toggle_menu}
                onClick={() => {
                  this.props.toggleHash("Help");
                }}
              >
                Help
              </span>{" "}
              menu.
            </p>
            <p>&nbsp;</p>
          </span>
        )}
      </div>
    );
  }
}
