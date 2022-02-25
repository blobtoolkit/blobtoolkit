import React from 'react'
import styles from './Fields.scss';

class FieldBoxHeaderButton extends React.Component {
  render(){
    return (
      <div className={styles.header_button} onClick={()=>{this.props.onAxisButtonClick(this.props.axis,this.props.fieldId)}}>
        {this.props.axis}
      </div>
    );
  }
}

export default FieldBoxHeaderButton
