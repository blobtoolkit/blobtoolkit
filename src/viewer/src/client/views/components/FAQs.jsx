import React from 'react'
import { connect } from 'react-redux'
import styles from './FAQ.scss'
import searchStyles from './Search.scss'
import SVGIcon from './SVGIcon'
import selectInverseIcon from './svg/invert.svg'
import xAxisIcon from './svg/xAxis.svg'
import categoryIcon from './svg/category.svg'

class FAQ extends React.Component {
  constructor(props) {
    super(props);
    this.state = {collapsed:true}
  }

  toggleState(){
    this.setState({collapsed:!this.state.collapsed})
  }

  render(){
    let css = this.state.collapsed ? styles.hidden : styles.answer
    let link
    if (this.props.link){
      link = <a href={this.props.url} target={this.props.link}>{this.props.link}</a>
    }
    let icon
    if (this.props.icon){
      icon = <SVGIcon sprite={this.props.icon} active={true}/>
    }
    return (
      <div className={styles.FAQ}>
        <div className={styles.question}
          onClick={()=>this.toggleState()}>
          <small>{this.state.collapsed ? '\u25B6' : '\u25BC'}</small> {this.props.question}
        </div>
        <div className={css}>
          {this.props.answer} {link} {icon}
        </div>
      </div>
    )
  }
}

export default class FAQs extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      query: ''
    }
  }

  handleInputChange(){
    if (this.search.value && this.search.value.length > 2){
      this.setState({
        query: this.search.value
      })
    }
    else {
      this.setState({query:undefined})
    }
  }

  render(){
    let questions = [
      {
        s: 'Datasets',
        q: 'What datasets are available?',
        a: 'Enter "all" in the search box in the "Datasets" menu to see all datasets.'
      },{
        s: 'Datasets',
        q: 'How do I find specific datasets?',
        a: 'Search for datasets by name, ID, phylum, etc. using the search box in the "Datasets" menu.'
      },{
        s: 'Datasets',
        q: 'How do I see the metadata for a dataset?',
        a: 'Metadata for an active dataset is summarised in the Detail view. Metadata for other datasets listed in the "Datasets" menu can be previewed by clicking on the Details button. Full metadata can be downloaded in JSON format from either of these views.'
      },{
        s: 'Datasets',
        q: 'How were the datasets analysed?',
        a: 'All datasets were analysed using the BlobToolKit INSDC pipeline. For more recent datasets the repository URL and commit hash for the version of the pipeline that was used are included in the metadata. The pipeline is available on github: ',
        u: 'https://github.com/blobtoolkit/insdc-pipeline',
        t: 'github.com/blobtoolkit/insdc-pipeline'
      },{
        s: 'Datasets',
        q: 'Can I load my own data?',
        a: 'This instance of the BlobToolKit Viewer is focussed on public datasets. To load your own data you will need to set up a local instance. Full instructions will be coming soon. For now visit the project github repository to get the code: ',
        u: 'https://github.com/blobtoolkit/viewer',
        t: 'github.com/blobtoolkit/viewer'
      },{
        s: 'Filters',
        q: 'How do I apply filters to a dataset?',
        a: 'Datasets can be filtered based on any category or variable using the "Filters" menu. To filter based on a category, click the tab above any bar of the preview histogram to hide/show records assigned to this category. To filter variables, use the sliders at either end of the preview histogram to change the maximum and minimum values or enter the numbers directly into the text boxes above the histogram.'
      },{
        s: 'Filters',
        q: 'How do I invert a filter?',
        a: 'Each filter can be inverted by clicking the invert icon: ',
        i: selectInverseIcon
      },{
        s: 'Filters',
        q: 'How do I filter based on a selection?',
        a: 'A selection can be used to filter a dataset by clicking the selection header in the "Filters" menu to activate it.',
      },{
        s: 'Selection',
        q: 'How do I select records?',
        a: 'Records can be selected by clicking bins in the Blobplot view or selecting individual records in the Table view.',
      },{
        s: 'Selection',
        q: 'How do I find out more about the records in a selection?',
        a: 'Once you have a selection, selected records are shown as a percentage of the total for each category in the "Summary" menu.',
      },{
        s: 'Plot',
        q: 'How do I change which variables are plotted?',
        a: 'To plot a variable on a given axis, use the "Filters" menu and click the relevant axis icon in the filter header: ',
        i: xAxisIcon
      },{
        s: 'Plot',
        q: 'How do I change which category is plotted?',
        a: 'Click the header of a category to activate it then click the category icon to use that category in the plots: ',
        i: categoryIcon
      },{
        s: 'Plot',
        q: 'What views are available?',
        a: 'Use the "Settings" menu to choose between Blobplot, Cumulative distibution, Detail, Snail plot and Table views. The Report view combines several views on a single page.',
      },{
        s: 'Plot',
        q: 'Can I display more than 9 categories?',
        a: 'No. The plots are limited to 9 categories to avoid having highly similar colours for different categories.',
      },{
        s: 'Plot',
        q: 'How do I see what is in the "other" category?',
        a: 'The category names that have been assigned to other are listed in the "Settings" menu.',
      },{
        s: 'Plot',
        q: 'Can I change the category order?',
        a: 'Almost. There is support for specifying a fixed category order (e.g. to plot a category individually that would otherwise be included in other), but this currently requires manually editing the URL.',
      },{
        s: 'Plot',
        q: 'How do I change the colours?',
        a: 'Two palettes are available in the "Settings" menu. To change colours click on one of the swatches in the user colour palette and select a new colour. To restore default colours, selct the default colour palette.',
      },{
        s: 'Plot',
        q: 'How do I change the bin size?',
        a: 'The size slider in the "Settings" menu can be used to make the bins larger or smaller.',
      },{
        s: 'Plot',
        q: 'What is the reducer function?',
        a: 'The reducer function (sum, maximum, minimum, count or mean) is applied to all z values for each category in a bin to obtain a single value to be plotted.',
      },{
        s: 'Plot',
        q: 'What is the scale function?',
        a: 'The relative sizes at which different z values are plotted is determined by the scale function (log10, linear or square-root).',
      },{
        s: 'Export',
        q: 'How do I save images?',
        a: 'Plots (and preview histograms from the "Filters" menu) can be saved as PNG or SVG files by clicking the buttons at the top right of the plot.',
      },{
        s: 'Export',
        q: 'How do I save a filtered view?',
        a: 'Filter parameters are added to the URL in the browser address bar so saving the URL provides a link back to a specific set of parameters. Selection-based filters cannot be stored in the URL so to export these you will need to use the "Lists" menu',
      },{
        s: 'Export',
        q: 'How do I export a filtered list of records?',
        a: 'Two lists are created by default in the "Lists" menu, one based on the current set of parameters (current) and a second unfiltered list (all) is available if the current list does not contain all records. The current list can be stored under a different name using the Create list button. To export a list of record names along with the parameter settings, click on the list then click Download JSON in the pop-up dialog box.',
      },{
        s: 'Export',
        q: 'Can I export raw data?',
        a: 'Data for all active variables can be exported in CSV format from the Table view by clicking the button at the top right of the table.',
      },{
        s: 'Import',
        q: 'How do I import a list?',
        a: 'Drag and drop a JSON file into the dotted upload area in the "Lists" menu. This function is designed to allow exported lists to be reloaded but lists could be created manually by copying the formatting of an exported list.',
      },{
        s: 'More',
        q: 'How do I cite this tool?',
        a: 'We are working on a manuscript describing this tool, the analysis pipeline and the datasets we have processed so far. Meanwhile, please cite the code using the Zenodo doi in the footer.',
      },{
        s: 'More',
        q: 'How do I report bugs/suggest features?',
        a: 'Please submit an issue at the project git repository: ',
        u: 'https://github.com/blobtoolkit/viewer',
        t: 'github.com/blobtoolkit/viewer'
      }
    ]
    let byCat = {Top:[]}
    let cat
    questions.forEach(q=>{
      let action = 'push'
      if (this.state.query && q.q.match(new RegExp('(^| )'+this.state.query,'i'))){
        cat = 'Top'
        action = 'unshift'
      }
      else if (this.state.query && q.a.match(new RegExp('(^| )'+this.state.query,'i'))){
        cat = 'Top'
      }
      else {
        cat = q.s
      }
      if (!byCat[cat]){
        byCat[cat] = []
      }
      if (action == 'push'){
        byCat[cat].push(q)
      }
      else {
        byCat[cat].unshift(q)
      }
    })
    let lists = []
    Object.keys(byCat).forEach(key=>{
      if (byCat[key].length > 0){
        if (key == 'Top'){
          lists.push(<span className={styles.result_count} key='span'>{byCat[key].length+' matching questions:'}</span>)
        }
        lists.push(
          <div className={styles.section} key={key}
            children={
              byCat[key].map((q,i)=>{
                return (
                  <FAQ key={i} question={q.q} answer={q.a} url={q.u} link={q.t} icon={q.i}/>
                )
              })
            }
          />
        )
        if (key == 'Top') lists.push(<hr key='hr'/>)
      }
    })
    let placeholder = "Search FAQs..."
    return (
      <form onSubmit={e=>e.preventDefault()}>
        <input
          placeholder={placeholder}
          ref={input => this.search = input}
          onChange={()=>this.handleInputChange()}
          className={searchStyles.search_box}
        />
        <div className={styles.outer}>
          {lists}
        </div>
      </form>

    )
  }
}
