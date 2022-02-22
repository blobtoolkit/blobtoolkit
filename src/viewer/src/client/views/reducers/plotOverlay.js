import { createSelector } from 'reselect'
import { getTransformFunction } from './plotParameters'
import { getTransformFunctionParams } from './plotParameters'
import { getSquareGrid } from './plotSquareBins'

export const getTransformLines = createSelector(
  getTransformFunctionParams,
  getTransformFunction,
  getSquareGrid,
  (params,transform,grid) => {
    let lines = []
    let origins = []
    let size = 900
    let limit = 900
    for (let j = 0; j <= size; j+= limit){
      let [x,y] = transform([0,j])
      let line = 'M0 '+(size-y)
      for (let i = 1; i <= size; i++){
        [x,y] = transform([i,j])
        line += ' L'+x+' '+(size-y)
      }
      lines.push(line)
    }
    // for (let i = 0; i <= size; i+= limit){
    //   let [x,y] = transform([i,0])
    //   let line = 'M'+x+' '+(size-y)
    //   for (let j = 1; j <= size; j++){
    //     [x,y] = transform([i,j])
    //     line += ' L'+x+' '+(size-y)
    //   }
    //   lines.push(line)
    // }
    let [x,y] = transform([0,limit - params.origin.y])
    let line = `M${x} ${limit - y}`;
    [x,y] = transform([limit,limit - params.origin.y])
    line += ` L${x} ${limit - y}`
    origins.push(line);
    [x,y] = transform([params.origin.x,0])
    line = `M${x} ${limit-y}`;
    [x,y] = transform([params.origin.x,limit])
    line += ` L${x} ${limit-y}`
    origins.push(line)
    return { lines, params, origins }
  }
)
