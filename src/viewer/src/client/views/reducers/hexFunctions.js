export const hexCorner = (center, radius, i) => {
  let angle_deg = 60 * i + 30
  let angle_rad = Math.PI / 180 * angle_deg
  return {
    cx:center.x + radius * Math.cos(angle_rad),
    cy:center.y + radius * Math.sin(angle_rad)
  }
}

export const cube_to_axial = cube => {
  let q = cube.x
  let r = cube.z
  return {q, r}
}

export const axial_to_cube = hex => {
  let x = hex.q
  let z = hex.r
  let y = -x-z
  return {x, y, z}
}

export const cube_to_oddr = cube => {
  let i = cube.x + (cube.z - (cube.z&1)) / 2
  let j = cube.z
  return {i, j}
}

export const oddr_to_cube = oddr => {
  let x = oddr.j - (oddr.i - (oddr.i&1)) / 2
  let z = oddr.i
  let y = -x-z
  return {x, y, z}
}

export const hex_round = hex => {
  return {q:Math.round(hex.q),r:Math.round(hex.r)}
}

export const pixel_to_hex = (x, y, radius) => {
  let q = (x * Math.sqrt(3)/3 - y / 3) / radius
  let r = y * 2/3 / radius
  return hex_round({q, r})
}

export const pixel_to_oddr = (x, y, radius) => {
  let q = Math.round((x * Math.sqrt(3)/3 - y / 3) / radius)
  let r = Math.round(y * 2/3 / radius)
  let i = q + (r - (r&1)) / 2
  let j = r
  return {i, j}
}

export const hex_to_oddr = hex => {
  let i = hex.q + (hex.r - (hex.r&1)) / 2
  let j = hex.r
  return {i, j}
}

export const drawPoints = (i,j,radius,res,scale = 1) => {
  let center
  let width = radius * Math.sqrt(3)
  let rowspacing = (3 / 2) * radius
  let corners = [0,1,2,3,4,5]
  if (j%2 == 0){
    center = {x:(i)*width,y:(j)*rowspacing}
  }
  else {
    center = {x:(i)*width+width/2,y:(j)*rowspacing}
  }
  if (center){
    radius = radius * scale
    let vertices = corners.map(i=>hexCorner(center,radius,i))
    let points = vertices.map(v=>v.cx+','+v.cy).join(' ')
    return points
  }
  return
}
