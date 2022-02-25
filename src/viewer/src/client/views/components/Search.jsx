import React, { Component } from "react";
import {
  fetchRepository,
  getAvailableDatasetIds,
} from "../reducers/repository";

import Suggestions from "./Suggestions";
import { apiUrl } from "../reducers/api";
import { connect } from "react-redux";
import styles from "./Search.scss";

// import { getSearhTerm,setSearchTerm } from '../reducers/search'


class SearchBox extends Component {
  constructor(props) {
    super(props);
    this.state = {
      query: "",
      results: [],
    };
  }

  getInfo() {
    fetch(apiUrl + "/search/autocomplete/" + this.state.query)
      .then(
        (response) => response.json(),
        (error) => console.log("An error occured.", error)
      )
      .then((json) => {
        this.setState({
          results: json.slice(0, 10),
        });
      });
    // let ret = this.props.onChangeText(this.state.query)
    // console.log(ret)
    //   .then(({ data }) => {
    //     this.setState({
    //       // results: data.data
    //       results: [{id:'test',name:data.data}]
    //     })
    //   })
  }

  handleInputChange() {
    this.setState(
      {
        query: this.search.value,
      },
      () => {
        if (this.state.query && this.state.query.length > 1) {
          // if (this.state.query.length % 2 === 0) {
          this.getInfo();
          // }
        } else {
          this.setState({ results: [] });
        }
      }
    );
  }

  handleKeyStroke(key, value) {
    if (key == "Enter") {
      if (this.state.results.length > 0) {
        this.props.onChooseTerm(this.state.results[0].term);
      }
    }
  }

  render() {
    let placeholder =
      this.props.datasetIds.length == 0
        ? "Search datasets to begin... (e.g. Nematoda)"
        : "Search for datasets...";
    return (
      <form onSubmit={(e) => e.preventDefault()}>
        <input
          placeholder={placeholder}
          ref={(input) => (this.search = input)}
          onChange={() => this.handleInputChange()}
          onKeyDown={(e) => this.handleKeyStroke(e.key, e.target.value)}
          className={styles.search_box}
        />
        <Suggestions
          results={this.state.results}
          onChooseTerm={this.props.onChooseTerm}
        />
      </form>
    );
  }
}

class Search extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = (state) => {
      return {
        datasetIds: getAvailableDatasetIds(state),
      };
    };
    this.mapDispatchToProps = (dispatch) => {
      return {
        onChangeText: (str) => {
          dispatch(fetchSearchResults(str));
        },
        onChooseTerm: (str) => {
          dispatch(fetchRepository(str));
        },
      };
    };
  }

  render() {
    const ConnectedSearchBox = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(SearchBox);
    return <ConnectedSearchBox {...this.props} />;
  }
}

export default Search;
