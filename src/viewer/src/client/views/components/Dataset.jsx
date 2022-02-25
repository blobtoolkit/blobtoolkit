import {
  getFieldHierarchy,
  getFieldsByParent,
  getTopLevelFields,
} from "../reducers/field";

import DatasetApplyFilters from "./DatasetApplyFilters";
import Field from "./Field";
import FieldSet from "./FieldSet";
import MainPlot from "./MainPlot";
import React from "react";
import Spinner from "./Spinner";
import { connect } from "react-redux";
import { getDatasetIsFetching } from "../reducers/repository";
import styles from "./Datasets.scss";

const mapStateToProps = (state) => {
  return {
    isFetching: getDatasetIsFetching(state),
    topLevelFields: getTopLevelFields(state),
    fields: getFieldHierarchy(state),
  };
};

const mapDispatchToProps = (dispatch) => {
  return {
    onDatasetClick: (id) => {},
  };
};

class Overview extends React.Component {
  mapFields(fields) {
    return fields.map((field) => {
      let jsx;
      if (field.hasRecords) {
        jsx = (
          <Field key={field.id} fieldId={field.id}>
            {field.id}
          </Field>
        );
      }
      if (field.children) {
        return (
          <FieldSet key={field.id + "_children"} title={field.id}>
            {jsx}
            {this.mapFields(field.children)}
          </FieldSet>
        );
      } else {
        return jsx;
      }
    });
  }

  render() {
    if (this.props.isFetching) {
      return <Spinner />;
    }
    let fields = this.mapFields(this.props.fields);
    return (
      <div className={styles.view_height}>
        <DatasetApplyFilters />
        {fields}
        <MainPlot />
      </div>
    );
  }
}

const Dataset = connect(mapStateToProps, mapDispatchToProps)(Overview);

export default Dataset;
