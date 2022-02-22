import JSONPretty from "react-json-pretty";
import Modal from "@material-ui/core/Modal";
import React from "react";
import styles from "./Layout.scss";

export class DatasetModal extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isShowingModal: Boolean(this.props.selected),
    };
  }
  handleClick() {
    this.setState({ isShowingModal: !this.state.isShowingModal });
  }
  handleClose() {
    this.setState({ isShowingModal: false });
    this.props.dismiss();
  }

  _downloadJSONFile(name, content) {
    var element = document.createElement("a");
    var file = new Blob([JSON.stringify(content)], { type: "text/plain" });
    element.href = URL.createObjectURL(file);
    element.download = (name || "file") + ".json";
    element.dispatchEvent(
      new MouseEvent(`click`, { bubbles: true, cancelable: true, view: window })
    );
  }

  render() {
    let meta = {};
    let full = {};
    let keys = ["id", "taxon", "assembly"];
    keys.forEach((k) => (meta[k] = this.props.meta[k]));
    if (this.props.meta.hasOwnProperty("revision")) {
      meta.revision = this.props.meta.revision;
    } else {
      meta.revision = 1;
    }
    Object.keys(this.props.meta).forEach((k) => {
      if (k != "fields" && k != "plot") {
        full[k] = this.props.meta[k];
      }
    });
    full.revision = meta.revision;
    return (
      <div
        style={{ position: "absolute", top: 0, right: 0, bottom: 0, left: 0 }}
        onClick={() => this.handleClick()}
      >
        <Modal
          open={this.state.isShowingModal}
          onClose={() => this.handleClose()}
        >
          <div className={styles.modalDialog}>
            <div className={styles.modal}>
              <a
                className={styles.button}
                onClick={() =>
                  this._downloadJSONFile(this.props.meta.id + ".meta", full)
                }
              >
                Download full metadata
              </a>
              <h2>{meta.name}</h2>
              <div className={styles.code_block}>
                <JSONPretty id="json-pretty" json={meta} />
              </div>
            </div>
          </div>
        </Modal>
      </div>
    );
  }
}

export default DatasetModal;
