import React from 'react'
import styles from './Layout.scss';
import Spinner from './Spinner'
import { format as d3Format } from 'd3-format'

const SelectedFraction = ({value,total}) => {
  return (
    <span>

    </span>
  )
}

class TabbedFraction extends React.Component {
  constructor(props) {
    super(props);
    this.state = {padTop: false}
  }

  handleDragEnter(e,title) {
    e.preventDefault()
    if (this.props.setTarget){
      if (!this.state.padTop){
        this.setState({padTop:true})
      }
      this.props.setTarget(title)
    }

  }
  handleDragLeave(e,title) {
    e.preventDefault()
    if (this.props.setTarget){
      this.setState({padTop:false})
      setTimeout(()=>{
        this.props.removeTarget(title)
      },100)
    }
  }

  // handleDragEnd(e) {
  //   e.preventDefault()
  //   if (this.state.padTop){
  //     this.setState({padTop:false})
  //     this.props.setTarget()
  //   }
  // }

  render() {
    let {value,total,color='rgba(0,0,0,0)',title='total',sel=false} = this.props
    let css = styles.colored_tab;
    value = Math.max(value,0)
    let portion = value / total
    let percent = d3Format(".1%")(portion)
    total = d3Format(",.0f")(total)
    value = d3Format(",.0f")(value)
    let style = this.state.padTop ? ({ paddingTop: this.props.padTop}) : {}
    return (
      <tr onDragEnter={(e)=>this.handleDragEnter(e,title)}
          onDragLeave={(e)=>this.handleDragLeave(e,title)}
          onDragOver={(e)=>{e.preventDefault()}}
          onDrop={(e)=>{e.preventDefault()}}
          draggable={true}
          onDragStart={(e)=>e.dataTransfer.setData('text','')}
          onDrag={(e)=>this.props.handleDragStart(e,title)}
          onDragEnd={(e)=>this.props.handleDragEnd(e)}>
        <td style={style}>
          <span
            className={css}
            style={{backgroundColor:color}}>
            &nbsp;
          </span>
        </td>
        <td className={styles.left_align} style={style}>
          {title}
        </td>
        <td style={style}>
          {sel ? (
            <span>
              <span className={styles.highlight}>
                {value}
              </span> /
            </span>
          ) : ' '}
        </td>
        <td style={style}>
          {total}
        </td>
        <td style={style}>
          {sel ? (
            <span className={styles.highlight}>
              {percent}
            </span>
          ) : ' '}
        </td>

      </tr>
    )
  }
}

class MenuSummary extends React.Component {
  constructor(props) {
    super(props);
    this.state = {source: false, target: false, title: false}
  }

  setTarget(title, target){
    if (title){
      this.setState({title, target})
    }
  }

  removeTarget(title, target){
    if (this.state.title == title){
      this.setState({title:false, target:false})
    }
  }

  handleDragStart(e, source){
    this.setState({source})
  }

  handleDragEnd(e){
    e.preventDefault()
    if (e.dataTransfer.dropEffect == 'move'
        && (this.state.target || this.state.target === 0)){
      this.changeOrder(this.state, this.props.bins)
    }
    this.setState({source: false, target: false, title: false})
  }

  changeOrder(state,bins){
    let values = {}
    let taxa = []
    let list
    bins = bins.slice(0)
    let counts = this.props.values.counts.binned
    let target = state.target
    let param = `${this.props.catAxis}--Order`
    let other = -1
    let current = this.props.parsed[param]
    if (current){
      let arr = current.split(',')
      let index = arr.indexOf(state.source)
      other = arr.indexOf('other')
      if (other > -1) other++
      if (index > -1){
        arr.splice(index,1)
      }
      if (arr.length >= state.target){
        let target = state.target
        while (target >= 0 && counts[target-1] == 0){
          target--
        }
        arr.splice(target,0,state.source)
        // if (index > -1){
        //   arr.splice(index,1)
        // }
        if (state.source == 'other'){
          other = target + 1
        }
        else {
          other++
        }
        list = arr.join(',')
      }
    }
    if (!list) {
      let index = bins.findIndex(o=>o.id==state.source)
      let target = index
      let skip = {}
      while (target >= 0 && counts[target-1] == 0){
        skip[target-1] = true
        target--
      }
      for (let i=0; i < state.target; i++){
        if (i != index && !skip[i]){
          taxa.push(bins[i].id)
        }
      }
      list = `${taxa.join(',')},${state.source}`
      if (state.source == 'other'){
        other = state.target + 1
      }
    }
    values[param] = list
    taxa = list.split(',')
    // let index = current.indexOf('other')
    // console.log(index)
    // let other
    // if (index != -1){
    //   other = index + 1
    //   values.otherLimit = other
    //   taxa = taxa.slice(0,index)
    //   values[param] = taxa.join(',')
    // }
    if (other > 0){
      if (other < 10){
        values.otherLimit = other
      }
      else {
        values.otherLimit = 10
      }
      // taxa = taxa.slice(0,other)
      // values[param] = taxa.join(',')
    }
    this.props.onChangeOrder(values,[])
    // this.props.changeQueryParams(Object.assign({[param]: list},this.props.parsed))
  }

  render() {
    let {values,zAxis,catAxis,bins,palette,other,fullSummary,parsed} = this.props
    let css = styles.long_menu_item
    let counts = []
    let reduced = []
    let listDiv
    let reset
    if (parsed[`${catAxis}--Order`]){
      reset = (<div className={styles.reset}
           onClick={()=>{this.props.onChangeOrder({},[`${catAxis}--Order`,'otherLimit'])}}
           >reset</div>)
    }
    if (bins){
      bins.forEach((bin,i) => {
        let title = bin.id
        let color = palette.colors[i]
        let value = values.counts.selBinned[i]
        let total = values.counts.binned[i]
        if (total){
          counts[i] = <TabbedFraction handleDrop={i > 0 && this.handleDrop}
                                      handleDragStart={(e,title)=>this.handleDragStart(e,title)}
                                      handleDragEnd={(e)=>this.handleDragEnd(e)}
                                      setTarget={(t)=>this.setTarget(t, i)}
                                      removeTarget={(t)=>this.removeTarget(t, i)}
                                      padTop={this.state.source ? '1.5em' : '0.5em'}
                                      index={i}
                                      key={i}
                                      {...{value,total,color,title,sel:values.reduced.sel}}/>
          value = values.reduced.selBinned[i]
          total = values.reduced.binned[i]
          reduced[i] = <TabbedFraction handleDrop={i > 0 && this.handleDrop}
                                       handleDragStart={(e,title)=>this.handleDragStart(e,title)}
                                       handleDragEnd={(e)=>this.handleDragEnd(e)}
                                       setTarget={(t)=>this.setTarget(t, i)}
                                       removeTarget={(t)=>this.removeTarget(t, i)}
                                       padTop={this.state.source ? '1.5em' : '0.5em'}
                                       index={i}
                                       key={i}
                                       {...{value,total,color,title,sel:values.counts.sel}}/>
        }
      })
      if (other && other.length){
        let list = []
        other.forEach((id,i) => {
          let css = styles.no_break
          list.push(
              <span key={i}>
                <span draggable={true}
                      className={css}
                      style={{cursor:'pointer'}}
                      onDragStart={(e)=>e.dataTransfer.setData('text','')}
                      onDrag={(e)=>this.handleDragStart(e,id)}
                      onDragEnd={(e)=>this.handleDragEnd(e)}>
                  {id}
                </span>, </span>
            )
        })
        listDiv = (<div className={css}>
          {reset}
          <h3>Other Taxa</h3>
          {list}
        </div>)
      }
    }
    return (
      <div>
        <div className={css}>
          {reset}
          <h3>{zAxis}</h3>
          <table className={styles.right_align}>
            <tbody>
              <TabbedFraction value={values.reduced.sel} total={values.reduced.all} sel={values.reduced.sel}/>
              {reduced}
            </tbody>
          </table>
        </div>
        <div className={css}>
          {reset}
          <h3>counts</h3>
          <table className={styles.right_align}>
            <tbody>
              <TabbedFraction value={values.counts.sel} total={values.counts.all} sel={values.counts.sel}/>
              {counts}
            </tbody>
          </table>
        </div>
        {listDiv}
      </div>
    )
  }
}

export default MenuSummary
