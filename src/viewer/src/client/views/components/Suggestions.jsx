import React from 'react'
import styles from './Search.scss';

const Suggestions = (props) => {
  const options = props.results.map((r,i) => {
    let terms
    if (r.names.length > 6){
      terms = r.names.slice(0,5).join(', ')+'... ['+r.names.length+' datasets]'
    }
    else {
      terms = r.names.join(', ')
    }
    return (
      <li
        onClick={()=>props.onChooseTerm(r.term)}
        key={r.term+'_'+r.field}>
        {r.term} <small>[{r.field}]</small> &ndash; {terms}
      </li>
    )
  })
  return <ul className={styles.autocomplete}>{options}</ul>
}

export default Suggestions
