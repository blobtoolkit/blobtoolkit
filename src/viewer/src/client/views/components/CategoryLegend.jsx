import React from 'react';
import styles from './Plot.scss'

export const CategoryLegend = ({categories,colors,otherColor}) => {
  let items = []
  if (categories){
    let w = 14
    let h = 14
    let gap = 110
    let offset = 900 - (w+gap)*Math.ceil(Object.keys(categories).length/2)
    let row = 0

    Object.keys(categories).forEach((key,i) => {
      let color = colors[key] || otherColor
      let title = categories[key] || 'extra'
      items.push(
          <g key={i} transform={'translate('+offset+','+(5+row*20)+')'}>
            <rect x={0} y={0} width={w} height={h} style={{fill:color,stroke:'black',strokeWidth:'1'}} />
            <text className={styles.horiz_legend} transform={'translate('+(w+3)+','+(h-2)+')'}>{title}</text>
          </g>
        )
        if (row){
          offset += w + gap
          row = 0
        }
        else {
          row = 1
        }
    })
  }
  return (
    <g>
      {items}
    </g>
  )
}

export default CategoryLegend
