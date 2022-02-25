import CategoryDistribution from "./CategoryDistribution";
import Modal from "@material-ui/core/Modal";
import React from "react";
import styles from "./Layout.scss";

export class RecordModal extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isShowingModal: Boolean(
        this.props.currentRecord || this.props.currentRecord === 0
      ),
    };
  }

  componentDidUpdate(nextProps) {
    if (this.props.currentRecord !== nextProps.currentRecord) {
      this.setState({
        isShowingModal: Boolean(
          this.props.currentRecord || this.props.currentRecord === 0
        ),
      });
    }
  }

  handleClick() {
    this.setState({ isShowingModal: true });
  }

  handleClose() {
    this.props.dismiss();
    this.setState({ isShowingModal: false });
  }

  render() {
    if (this.state.isShowingModal) {
      return (
        <div className={styles.modal}>
          <Modal
            open={this.state.isShowingModal}
            onClose={() => this.handleClose()}
          >
            <div className={styles.modalDialog}>
              <CategoryDistribution />
            </div>
          </Modal>
        </div>
      );
    }
    return <div className={styles.modal} />;
  }
}

export default RecordModal;
