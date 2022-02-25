import * as d3 from "d3";

import { Link } from "react-router-dom";
import React from "react";
import ReactDOM from "react-dom";
import Resizable from "re-resizable";
import styles from "./Filters.scss";
import testSprite from "./svg/test.svg";
const apiUrl = "/api/v1";

class AvailableFiltersBox extends React.Component {
  render() {
    var children = [];
    this.props.fields.forEach((field) => {
      children.push(
        <FilterBox
          datasetId={this.props.datasetId}
          filter={field.filters["default"]}
          key={field._id}
        />
      );
    });
    return <div>{children}</div>;
  }
}

class FilterBox extends React.Component {
  constructor(props) {
    super(props);
    let x;
    if (props.filter.scale) {
      x = props.filter.scale;
      //if (this.props.filter.id == 'gc'){ x = d3.scaleLinear() }
      x.domain(this.props.filter._limits).range([0, 400]); // TODO: remove magic numbers
    }
    this.state = {
      range: props.filter._limits,
      xScale: x,
    };
  }
  updateRange(value, index) {
    var newRange = this.state.range.slice(0);
    newRange[index] = value * 1;
    this.setState({
      range: newRange,
    });
  }
  xScale() {
    var newRange = this.state.range.slice(0);
    newRange[index] = value * 1;
    this.setState({
      range: newRange,
    });
  }
  render() {
    let headerComponents = [];
    let handleComponents = [];
    handleComponents.push(
      <FilterHandle
        key="right"
        handlePosition="right"
        range={this.state.range}
        update={this.updateRange.bind(this)}
        xScale={this.state.xScale}
      />
    );
    handleComponents.push(
      <FilterHandle
        key="left"
        handlePosition="left"
        range={this.state.range}
        update={this.updateRange.bind(this)}
        xScale={this.state.xScale}
      />
    );
    let preview;
    let handles;
    if (this.props.filter.field._preload) {
      headerComponents.push(
        <FilterRange
          key="range"
          range={this.state.range}
          update={this.updateRange.bind(this)}
        />
      );
      preview = (
        <FilterDataPreview
          datasetId={this.props.datasetId}
          filter={this.props.filter}
          xScale={this.state.xScale}
        />
      );
      handles = <FilterHandles children={handleComponents} />;
    }
    headerComponents.push(
      <FilterControls
        key="controls"
        datasetId={this.props.datasetId}
        filter={this.props.filter}
      />
    );

    return (
      <div id={this.props.filter.id} className={styles.outer}>
        <FilterHeader
          name={this.props.filter.field._name}
          children={headerComponents}
        />
        <div className={styles.main}>
          {preview}
          {handles}
        </div>
      </div>
    );
  }
}

const withLogProps = (MyComponent) =>
  class extends MyComponent {
    render() {
      // Wraps the input component in a container, without mutating it. Good!
      return <MyComponent {...this.props} />;
    }
  };
const EnhancedComponent = withLogProps(FilterBox);

// class FilterHeader extends React.Component {
//   render() {
//     return (
//       <div className={styles.header}>
//         <h1>{this.props.filter.name}</h1>
//         <FilterRange datasetId={this.props.datasetId} filter={this.props.filter}/>
//         <FilterControls datasetId={this.props.datasetId} filter={this.props.filter}/>
//       </div>
//     )
//   }
// }

const FilterHeader = ({ name, children }) => (
  <div className={styles.header}>
    <h1>{name}</h1>
    {children}
  </div>
);

class FilterRange extends React.Component {
  render() {
    return (
      <div className={styles.range}>
        <input
          type="text"
          className={styles.range_input}
          value={this.props.range[0]}
          onChange={(e) => {
            this.props.update(e.target.value, 0);
          }}
        />
        &nbsp;:&nbsp;
        <input
          type="text"
          className={styles.range_input}
          value={this.props.range[1]}
          onChange={(e) => {
            this.props.update(e.target.value, 1);
          }}
        />
      </div>
    );
  }
}
// class FilterRange extends React.Component {
//   render() {
//     return (
//       <div className={styles.range}>
//         <input type='text' className={styles.range_input} defaultValue={this.props.filter.limits[0]}/>
//         &nbsp;:&nbsp;
//         <input type='text' className={styles.range_input} defaultValue={this.props.filter.limits[1]}/>
//       </div>
//     )
//   }
// }

class FilterControls extends React.Component {
  render() {
    return (
      <div className={styles.controls}>
        <FilterSwitch
          datasetId={this.props.datasetId}
          filter={this.props.filter}
        />
      </div>
    );
  }
}

class FilterSwitch extends React.Component {
  render() {
    return (
      <svg viewBox={testSprite.viewBox} className={styles.icon}>
        <use xlinkHref={"#" + testSprite.id} />
      </svg>
    );
  }
}

class FilterDataPreview extends React.Component {
  constructor(props) {
    super(props);
    this.state = { loading: true };
  }

  componentDidMount() {
    d3.json(
      apiUrl +
        "/field/" +
        this.props.datasetId +
        "/" +
        this.props.filter.field._id,
      (error, data) => {
        if (error) {
          console.error(error);
        } else {
          this.state.data = data.values;
          this.drawChart();
        }
      }
    );
  }

  shouldComponentUpdate() {
    return false;
  }

  componentWillUnmount() {
    // ReactDOM.unmountComponentAtNode(this.tooltipTarget);
  }

  drawChart() {
    /*
      D3 code to create our visualization by appending onto this.svg
    */

    var svg = d3.select(this.svg);

    var height = this.svg.clientHeight;
    var width = this.svg.clientWidth;

    var data = this.state.data; //d3.range(1000).map(d3.randomBates(10));
    var g = svg.append("g");

    var x = this.props.xScale;
    var thresh = Array.from(Array(24).keys()).map((n) => {
      return x.invert(((n + 1) * width) / 25);
    });
    var bins = d3.histogram().domain(x.domain()).thresholds(thresh)(data);
    var y = d3
      .scaleLinear()
      .domain([
        0,
        d3.max(bins, function (d) {
          return d.length;
        }),
      ])
      .range([height, 0]);

    var bar = g
      .selectAll("." + styles.bar)
      .data(bins)
      .enter()
      .append("g")
      .attr("class", styles.bar)
      .attr("transform", function (d) {
        return "translate(" + x(d.x0) + "," + y(d.length) + ")";
      });

    bar
      .append("rect")
      .attr("x", 1)
      .attr("width", function (d) {
        return x(d.x1) - x(d.x0) - 1;
      })
      .attr("height", function (d) {
        return height - y(d.length);
      });
    // // At some point we render a child, say a tooltip
    // const tooltipData = ...
    // this.renderTooltip([50, 100], tooltipData);
  }

  render() {
    return (
      <div className={styles.data_preview_container}>
        {
          //<div ref={(elem) => { this.tooltipTarget = elem; }} />
        }
        <svg
          ref={(elem) => {
            this.svg = elem;
          }}
        ></svg>
      </div>
    );
  }

  // renderTooltip(coordinates, tooltipData) {
  //   const tooltipComponent = (
  //     <MyTooltipComponent
  //       coordinates={coordinates}
  //       data={tooltipData} />
  //   );
  //
  //   ReactDOM.render(tooltipComponent, this.tooltipTarget);
  // }
}

class FilterHandles extends React.Component {
  render() {
    return (
      <div className={styles.resizables_container}>{this.props.children}</div>
    );
  }
}

class FilterHandle extends React.Component {
  bound() {
    return this.props.handlePosition == "right" ? 1 : 0;
  }
  boundPercent() {
    return this.props.xScale(this.props.range[this.bound()]) / 4;
  }
  render() {
    return (
      <Resizable
        className={styles.resizable}
        width={this.boundPercent() + "%"}
        bounds={"parent"}
        height={"100%"}
        handlerClasses={{
          right: `${styles.handle} ${styles[this.props.handlePosition]}`,
        }}
        enable={{
          top: false,
          right: true,
          bottom: false,
          left: false,
          topRight: false,
          bottomRight: false,
          bottomLeft: false,
          topLeft: false,
        }}
        onResize={(e, direction, ref, delta) => {}}
        onResizeStop={(e, direction, ref, delta) => {
          this.props.update(
            this.props.xScale.invert(ref.clientWidth),
            this.bound()
          );
        }}
      />
    );
  }
}

export default AvailableFiltersBox;

export { AvailableFiltersBox, FilterBox };
