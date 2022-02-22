import React from 'react'
import styles from './Fields.scss';
import plotStyles from './Plot.scss';
import FieldBoxHeader from './FieldBoxHeader'
import FieldRawDataPreview from './FieldRawDataPreview'
import FilterPreview from './FilterPreview'
import Filter from './Filter'
import { scaleLinear as d3scaleLinear } from 'd3-scale';
import { scaleSqrt as d3scaleSqrt } from 'd3-scale';
import PreviewPlotBoundary from './PreviewPlotBoundary'

class FieldBox extends React.Component {
  constructor(props) {
    super(props);
    this.state = {yScale:d3scaleSqrt().domain([0, 100]).range([100, 0])}
  }
  toggleState(key){
    let obj = {id:this.props.fieldId};
    obj[key] = this.props.hasOwnProperty(key) ? !this.props[key] : true;
    this.props.toggleActive(obj);
    this.props.showData(this.props.fieldId);
  }
  render(){
    let outer_css = styles.outer
    let dimensions = {width:400,height:400/3.75}
    let viewbox = '0 0 '+(dimensions.width+150)+' '+(dimensions.height+80)
    let fieldContent
    let filterContent
    let filterPreview
    let css = styles.main
    let svgCss = styles.preview_svg
    if (this.props.type == 'category'){
      css += ' '+styles.cat
      svgCss += ' '+styles.cat
    }
    if (this.props.active && ! this.props.isStatic && !this.props.fieldId.match('selection')){
      outer_css += ' '+(this.props.type == 'category' ? styles.cat_expanded : styles.expanded);
      filterContent = <Filter {...this.props} />
      filterPreview = <FilterPreview {...this.props} dimensions={dimensions}/>
      fieldContent = (
        <div className={css}>
          <div className={svgCss}>
            <svg id={this.props.fieldId+'.preview'}
                viewBox={viewbox}
                preserveAspectRatio='xMinYMax meet'
                >
              <g transform='translate(75,5)'>
                <FieldRawDataPreview {...this.props} updateYScale={(obj)=>{this.setState(obj)}}  dimensions={dimensions}/>
                {filterPreview}
                <PreviewPlotBoundary xLabel={this.props.fieldId} yScale={this.state.yScale} xScale={this.props.xScale}  dimensions={dimensions}/>
              </g>
            </svg>
          </div>
        </div>
      )
    }

    let filterType = false;
    return (
      <div id={this.props.fieldId} className={outer_css} onClick={()=>{}}>
        <FieldBoxHeader {...this.props} onHeaderClick={()=>{this.toggleState('active')}}
         onAxisButtonClick={(axis,id)=>{this.props.setAxes(axis,id)}}
         onSumButtonClick={(id)=>{this.props.sumField(id)}}
         onInvertButtonClick={(id)=>{this.props.toggleInvert(id)}}/>
        {fieldContent}
        {filterContent}
      </div>
    );
  }
}

export default FieldBox
