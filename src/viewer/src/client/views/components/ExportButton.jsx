import React from 'react'
import { connect } from 'react-redux'
import styles from './Plot.scss'
const saveSvgAsPng = require('save-svg-as-png/lib/saveSvgAsPng.js')
import { getPngResolution } from '../reducers/plotParameters'


const _downloadJSONFile = (name,content) => {
  var element = document.createElement('a');
  var file = new Blob([JSON.stringify(content)], {type: 'text/plain'});
  element.href = URL.createObjectURL(file);
  element.download = (name || 'file') + '.json';
  element.dispatchEvent(new MouseEvent(`click`, {bubbles: true, cancelable: true, view: window}))
}

const _downloadTextFile = (name,content,format='text') => {
  var element = document.createElement('a');
  var file = new Blob([content+"\n"], {type: 'text/'+format});
  element.href = URL.createObjectURL(file);
  element.download = (name || 'file') + '.'+format;
  element.dispatchEvent(new MouseEvent(`click`, {bubbles: true, cancelable: true, view: window}))
}

const arrayBufferToBase64 = buffer => {
  let binary = '';
  let bytes = [].slice.call(new Uint8Array(buffer));
  bytes.forEach((b) => binary += String.fromCharCode(b));
  return window.btoa(binary);
};

const fetchImage = async url => {
  let response = await fetch(url)
  if (!response.ok){
    // console.log(response)
    return
  }
  let buffer = await response.arrayBuffer()
  return buffer;
}

const _downloadImgFile = async (name, url, format) => {
  let element = document.createElement('a');
  let content = await fetchImage(url)
  let file = new Blob([content], {type: 'image/'+(format == 'svg' ? 'svg+xml' : format)});
  element.href = URL.createObjectURL(file);
  element.download = (name || 'file') + '.'+format;
  element.dispatchEvent(new MouseEvent(`click`, {bubbles: true, cancelable: true, view: window}))

}

const screenScale = () => {
  let query = "(-webkit-min-device-pixel-ratio: 2), (min-device-pixel-ratio: 2), (min-resolution: 192dpi)";
  if (window.matchMedia(query).matches) {
    return 2
  }
  else {
    return 1
  }
}

const Button = ({view,element,data,url,format,prefix,scale=1,func,size=1000,res=1000}) => {
  if (!func){
    func = () => {}
    if (format == 'svg'){
      if (element){
        func = ()=>(saveSvgAsPng.saveSvg(document.getElementById(element),prefix+'.svg'))
      }
      else {
        func = ()=>_downloadImgFile(prefix,url,'svg')
      }
    }
    else if (format == 'png'){
      if (element){
        let scale = res / size / screenScale()
        func = ()=>(saveSvgAsPng.saveSvgAsPng(document.getElementById(element),prefix+'.png',{backgroundColor:'white',scale}))
      }
      else {
        func = ()=>_downloadImgFile(prefix,url,'png')
      }
    }
    else if (format == 'json'){
      func = ()=>_downloadJSONFile(prefix,data)
    }
    else {
      func = ()=>_downloadTextFile(prefix,data,format)
    }
  }
  return (
    <a id={`${view}_save_${format}`} className={styles.save_svg} onClick={func}>&#8681;{format}</a>
  )
}

export class ExportButton extends React.Component {
  constructor(props) {
    super(props);
    this.mapStateToProps = state => {
      return {
        res: getPngResolution(state)
      }
    }
  }

  render(){
    const ConnectedExportButton = connect(
      this.mapStateToProps
    )(Button)
    return (
      <ConnectedExportButton {...this.props}/>
    )
  }
}

export default ExportButton
