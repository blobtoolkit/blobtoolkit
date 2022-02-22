import React from 'react'
import styles from './Layout.scss';
import Dropzone from 'react-dropzone'

const ListUpload = ({active,onDrop}) => {
  return (
    <div className={styles.dropzone_menu_item}>
      <Dropzone
        accept="application/json"
        onDrop={onDrop}
        className={styles.dropzone}
        activeClassName={styles.dropzone_active}
        acceptClassName={styles.dropzone_accept}
        rejectClassName={styles.dropzone_reject}
      >
        {({ isDragActive, isDragReject, acceptedFiles, rejectedFiles }) => {
          if (isDragActive) {
            return "Drop file here to upload";
          }
          if (isDragReject) {
            return "This file is not authorized";
          }
          return acceptedFiles.length || rejectedFiles.length
            ? `Accepted ${acceptedFiles.length}, rejected ${rejectedFiles.length} files`
            : "Drag & drop or click to upload files (JSON format only)";
        }}
      </Dropzone>
    </div>
  )
}

export default ListUpload
