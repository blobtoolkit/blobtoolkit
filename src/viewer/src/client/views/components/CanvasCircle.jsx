import React from 'react';
import ReactDOM from 'react-dom';
import {Layer, Circle, Stage, Group} from 'react-konva';



class CanvasCircle extends React.Component {

    render() {
        return (
            <Circle
                x={this.props.x} y={this.props.y} radius={this.props.z}
                fill={this.props.color}
                stroke='#999'
                strokeWidth={0.25}
            />
        );
    }
}

export default CanvasCircle
