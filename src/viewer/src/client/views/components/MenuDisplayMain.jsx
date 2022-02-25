import {
  chooseCircleLimit,
  chooseCurveOrigin,
  chooseErrorBars,
  chooseLargeFonts,
  chooseMaxCount,
  chooseMaxSpan,
  choosePlotGraphics,
  choosePlotResolution,
  choosePlotScale,
  choosePlotShape,
  choosePngResolution,
  chooseSVGThreshold,
  chooseScaleTo,
  chooseShowTotal,
  chooseSideMax,
  chooseSnailOrigin,
  chooseWindowSize,
  chooseZReducer,
  chooseZScale,
  getCircleLimit,
  getCurveOrigin,
  getErrorBars,
  getLargeFonts,
  getMaxCount,
  getMaxSpan,
  getPlotGraphics,
  getPlotResolution,
  getPlotScale,
  getPlotShape,
  getPngResolution,
  getSVGThreshold,
  getScaleTo,
  getShowTotal,
  getSideMax,
  getSnailOrigin,
  getTransformFunctionParams,
  getWindowSize,
  getZReducer,
  getZScale,
  setPlotResolution,
  setTransformFunction,
} from "../reducers/plotParameters";
import {
  chooseNohitThreshold,
  chooseStaticThreshold,
  getNohitThreshold,
  getStaticThreshold,
} from "../reducers/repository";
import {
  chooseView,
  getDatasetID,
  getParsedQueryString,
  getStatic,
  getView,
  setQueryString,
  toggleStatic,
} from "../reducers/location";
import { getBuscoSets, getRecordCount } from "../reducers/summary";
import { getSelectedDatasetMeta, getStaticFields } from "../reducers/dataset";

import MenuDataset from "./MenuDataset";
import MenuDisplayKite from "./MenuDisplayKite";
import MenuDisplaySet from "./MenuDisplaySet";
import MenuDisplaySimple from "./MenuDisplaySimple";
import MenuItem from "./MenuItem";
import MenuSwitchView from "./MenuSwitchView";
import NumericInput from "./NumericInput";
import Palettes from "../containers/Palettes";
import React from "react";
import ReferenceAssembliesList from "./ReferenceAssembliesList";
import SVGIcon from "./SVGIcon";
import TextIcon from "./TextIcon";
import Timeout from "./Timeout";
import ToolTips from "./ToolTips";
import adjustIcon from "./svg/adjust.svg";
import centerIcon from "./svg/center.svg";
import circleIcon from "./svg/circleShape.svg";
import { connect } from "react-redux";
import countIcon from "./svg/count.svg";
import { format as d3Format } from "d3-format";
import { format as d3format } from "d3-format";
import fontSizeIcon from "./svg/fontSize.svg";
//import styles from './Plot.scss'
import { getAxisTitle } from "../reducers/plotData";
import { getFields } from "../reducers/field";
import { getMainPlotData } from "../reducers/plotData";
import hexIcon from "./svg/hexShape.svg";
import invertIcon from "./svg/invert.svg";
import kiteIcon from "./svg/kiteShape.svg";
import lineIcon from "./svg/lineShape.svg";
import linearIcon from "./svg/linear.svg";
import logIcon from "./svg/log.svg";
import maxIcon from "./svg/max.svg";
import meanIcon from "./svg/mean.svg";
import minIcon from "./svg/min.svg";
import { queryToStore } from "../querySync";
import scaleIcon from "./svg/scale.svg";
import sqrtIcon from "./svg/sqrt.svg";
import squareIcon from "./svg/squareShape.svg";
import styles from "./Layout.scss";
import sumIcon from "./svg/sum.svg";
import svgIcon from "./svg/svg.svg";
import xIcon from "./svg/xLetter.svg";
import yIcon from "./svg/yLetter.svg";
import zeroIcon from "./svg/zero.svg";

let pct = (n) => d3Format(".0%")(n);
let sci = (n) => d3Format(".1s")(n);

class DisplayMenu extends React.Component {
  render() {
    let {
      datasetId,
      title,
      view,
      busco,
      records,
      isStatic,
      hasStatic,
      meta,
      shape,
      onSelectShape,
      resolution,
      onChangeResolution,
      graphics,
      onChangeGraphics,
      threshold,
      onChangeThreshold,
      circleLimit,
      onChangeCircleLimit,
      reducer,
      onSelectReducer,
      scale,
      onSelectScale,
      curveOrigin,
      onSelectCurveOrigin,
      scaleTo,
      onSelectScaleTo,
      snailOrigin,
      onSelectSnailOrigin,
      transform,
      onChangeTransform,
      fields,
      plotScale,
      onChangePlotScale,
      onSelectView,
      onToggleStatic,
      pngResolution,
      onChangePngResolution,
      staticThreshold,
      onChangeStaticThreshold,
      nohitThreshold,
      onChangeNohitThreshold,
      sideMax,
      onChangeSideMax,
      data = {},
      onChangeAxisRange,
      parsed,
      changeQueryParams,
      showTotal,
      onChangeShowTotal,
      maxSpan,
      onChangeMaxSpan,
      maxCount,
      onChangeMaxCount,
      // adjustCoverage, onToggleAdjustCoverage,
      largeFonts,
      onToggleLargeFonts,
      errorBars,
      onChangeErrorBars,
      windowSize,
      onChangeWindowSize,
    } = this.props;
    let context;
    view = view || "blob";
    let blob;
    let xMeta = data.meta ? data.meta.x.meta : { limit: [0, 1], range: [0, 1] };
    let yMeta = data.meta ? data.meta.y.meta : { limit: [0, 1], range: [0, 1] };
    const resetLimits = (meta) => {
      let values = {};
      if (meta.limit[0] != meta.range[0]) {
        values[meta.id + "--LimitMin"] = meta.range[0];
        values[meta.id + "--Min"] = meta.range[0];
      }
      if (meta.limit[1] != meta.range[1]) {
        values[meta.id + "--LimitMax"] = meta.range[1];
        values[meta.id + "--Max"] = meta.range[1];
      }
      let remove = [
        meta.id + "--LimitMin",
        meta.id + "--Min",
        meta.id + "--LimitMax",
        meta.id + "--Max",
      ];
      remove.forEach((key) => {
        delete parsed[key];
      });
      onChangeAxisRange(values, remove);
      changeQueryParams(parsed);
    };
    let displayTotal = (
      <MenuDisplaySimple name="display total">
        <div className={styles.full_height}>
          <TextIcon
            title={"show"}
            active={showTotal}
            onIconClick={() => onChangeShowTotal("true")}
          />
          <TextIcon
            title={"hide"}
            active={!showTotal}
            onIconClick={() => onChangeShowTotal("false")}
          />
        </div>
      </MenuDisplaySimple>
    );
    if (isStatic) {
      blob = (
        <MenuDisplaySimple name="shape">
          <SVGIcon
            sprite={squareIcon}
            active={shape == "square"}
            onIconClick={() => onSelectShape("square")}
          />
          <SVGIcon
            sprite={hexIcon}
            active={shape == "hex"}
            onIconClick={() => onSelectShape("hex")}
          />
          <SVGIcon
            sprite={circleIcon}
            active={shape == "circle"}
            onIconClick={() => onSelectShape("circle")}
          />
        </MenuDisplaySimple>
      );
    } else {
      let windowSizes, errorBarOptions;
      if (meta.settings && meta.settings.stats_windows) {
        windowSizes = [];
        meta.settings.stats_windows.forEach((value) => {
          if (value < 1) {
            windowSizes.push({ label: pct(value), value });
          } else if (value > 1) {
            windowSizes.push({ label: sci(value), value });
          }
        });
        errorBarOptions = [
          { label: "0", value: 0 },
          { label: "SD", value: "sd" },
          { label: "SE", value: "se" },
          { label: "CI", value: "ci" },
        ];
      }
      blob = (
        <span>
          <MenuDisplaySimple name="shape">
            <SVGIcon
              sprite={squareIcon}
              active={shape == "square"}
              onIconClick={() => onSelectShape("square")}
            />
            <SVGIcon
              sprite={hexIcon}
              active={shape == "hex"}
              onIconClick={() => onSelectShape("hex")}
            />
            <SVGIcon
              sprite={circleIcon}
              active={shape == "circle"}
              onIconClick={() => onSelectShape("circle")}
            />
            {fields.gc_windows && records <= threshold && (
              <SVGIcon
                sprite={lineIcon}
                active={shape == "lines"}
                onIconClick={() => onSelectShape("lines")}
              />
            )}
            <SVGIcon
              sprite={kiteIcon}
              active={shape == "kite"}
              onIconClick={() => onSelectShape("kite")}
            />
          </MenuDisplaySimple>
          {shape == "lines" && windowSizes && (
            <MenuDisplaySimple name="window size">
              {windowSizes.map((obj) => (
                <TextIcon
                  key={obj.value}
                  title={obj.label}
                  active={obj.value == windowSize}
                  onIconClick={() => onChangeWindowSize(obj.value)}
                />
              ))}
            </MenuDisplaySimple>
          )}
          {shape == "lines" && errorBarOptions && (
            <MenuDisplaySimple name="error bars">
              {errorBarOptions.map((obj) => (
                <TextIcon
                  key={obj.value}
                  title={obj.label}
                  active={obj.value == errorBars}
                  onIconClick={() => onChangeErrorBars(obj.value)}
                />
              ))}
            </MenuDisplaySimple>
          )}
          {shape == "kite" && <MenuDisplayKite />}
          <MenuDisplaySimple name={"resolution [ " + resolution + " ]"}>
            <div className={styles.full_height}>
              <span className={styles.middle}>50</span>
              <input
                onChange={(e) => onChangeResolution(e.target.value)}
                type="range"
                value={resolution}
                min="5"
                max="50"
                step="1"
                className={styles.flip_horiz + " " + styles.middle}
                data-tip
                data-for="size-slider"
              />
              <span className={styles.middle}>5</span>
            </div>
          </MenuDisplaySimple>
          <MenuDisplaySimple name="reducer function">
            <SVGIcon
              sprite={sumIcon}
              active={reducer.id == "sum"}
              onIconClick={() => onSelectReducer("sum")}
            />
            <SVGIcon
              sprite={maxIcon}
              active={reducer.id == "max"}
              onIconClick={() => onSelectReducer("max")}
            />
            <SVGIcon
              sprite={minIcon}
              active={reducer.id == "min"}
              onIconClick={() => onSelectReducer("min")}
            />
            <SVGIcon
              sprite={countIcon}
              active={reducer.id == "count"}
              onIconClick={() => onSelectReducer("count")}
            />
            <SVGIcon
              sprite={meanIcon}
              active={reducer.id == "mean"}
              onIconClick={() => onSelectReducer("mean")}
            />
          </MenuDisplaySimple>
          <MenuDisplaySimple name="scale function">
            <SVGIcon
              sprite={logIcon}
              active={scale == "scaleLog"}
              onIconClick={() => onSelectScale("scaleLog")}
            />
            <SVGIcon
              sprite={linearIcon}
              active={scale == "scaleLinear"}
              onIconClick={() => onSelectScale("scaleLinear")}
            />
            <SVGIcon
              sprite={sqrtIcon}
              active={scale == "scaleSqrt"}
              onIconClick={() => onSelectScale("scaleSqrt")}
            />
          </MenuDisplaySimple>
          {/* <MenuDisplaySimple name='adjust coverage'>
            <SVGIcon sprite={adjustIcon} active={adjustCoverage} onIconClick={()=>onToggleAdjustCoverage(!adjustCoverage)}/>
          </MenuDisplaySimple> */}
          <MenuDisplaySimple
            name={"scale factor [ " + d3Format(",.1f")(plotScale) + " ]"}
          >
            <div className={styles.full_height}>
              <span className={styles.middle}>0.1</span>
              <input
                onChange={(e) => onChangePlotScale(e.target.value)}
                type="range"
                value={plotScale}
                min="0.1"
                max="2"
                step="0.1"
                className={styles.middle}
                data-tip
                data-for="scale-slider"
              />
              <span className={styles.middle}>2.0</span>
            </div>
          </MenuDisplaySimple>
          <MenuDisplaySimple name="x-axis range">
            <div className={styles.full_height}>
              {(xMeta.limit[0] != xMeta.range[0] ||
                xMeta.limit[1] != xMeta.range[1]) && (
                <span
                  className={styles.reset}
                  onClick={() => resetLimits(xMeta)}
                >
                  reset
                </span>
              )}
              <NumericInput
                initialValue={
                  isNaN(xMeta.limit[0]) ? xMeta.range[0] : xMeta.limit[0]
                }
                onChange={(value) =>
                  onChangeAxisRange({
                    [xMeta.id + "--LimitMin"]: value,
                    [xMeta.id + "--Min"]: value,
                  })
                }
              />
              <NumericInput
                initialValue={
                  isNaN(xMeta.limit[1]) ? xMeta.range[1] : xMeta.limit[1]
                }
                onChange={(value) =>
                  onChangeAxisRange({
                    [xMeta.id + "--LimitMax"]: value,
                    [xMeta.id + "--Max"]: value,
                  })
                }
              />
            </div>
          </MenuDisplaySimple>
          {yMeta && (
            <MenuDisplaySimple name="y-axis range">
              <div className={styles.full_height}>
                {yMeta.limit &&
                  yMeta.range &&
                  (yMeta.limit[0] != yMeta.range[0] ||
                    yMeta.limit[1] != yMeta.range[1]) && (
                    <span
                      className={styles.reset}
                      onClick={() => resetLimits(yMeta)}
                    >
                      reset
                    </span>
                  )}
                <NumericInput
                  initialValue={
                    isNaN(yMeta.limit[0]) ? yMeta.range[0] : yMeta.limit[0]
                  }
                  onChange={(value) =>
                    onChangeAxisRange({
                      [yMeta.id + "--LimitMin"]: value,
                      [yMeta.id + "--Min"]: value,
                    })
                  }
                />
                <NumericInput
                  initialValue={
                    isNaN(yMeta.limit[1]) ? yMeta.range[1] : yMeta.limit[1]
                  }
                  onChange={(value) =>
                    onChangeAxisRange({
                      [yMeta.id + "--LimitMax"]: value,
                      [yMeta.id + "--Max"]: value,
                    })
                  }
                />
              </div>
            </MenuDisplaySimple>
          )}
          <ReferenceAssembliesList />
          {shape == "circle" && (
            <MenuDisplaySimple name="plot graphics">
              <div className={styles.full_height}>
                <label className={styles.middle} htmlFor="#svgThreshold">
                  Threshold
                </label>
                <NumericInput
                  initialValue={threshold}
                  onChange={onChangeThreshold}
                />
                <span className={styles.middle}>or</span>
              </div>
              <SVGIcon
                sprite={svgIcon}
                active={graphics == "svg"}
                onIconClick={() =>
                  onChangeGraphics(graphics == "svg" ? "canvas" : "svg")
                }
              />
            </MenuDisplaySimple>
          )}
          {shape == "circle" && (
            <MenuDisplaySimple name="circle limit">
              <div className={styles.full_height}>
                <NumericInput
                  initialValue={circleLimit}
                  onChange={onChangeCircleLimit}
                />
              </div>
            </MenuDisplaySimple>
          )}
          <MenuDisplaySimple name="histogram maximum">
            <div className={styles.full_height}>
              <NumericInput initialValue={sideMax} onChange={onChangeSideMax} />
            </div>
          </MenuDisplaySimple>
          {displayTotal}
        </span>
      );
    }
    let cumulative = (
      <span>
        <MenuDisplaySimple name="curve origin">
          <SVGIcon
            sprite={zeroIcon}
            active={curveOrigin == "0"}
            onIconClick={() => onSelectCurveOrigin("0")}
          />
          <SVGIcon
            sprite={xIcon}
            active={curveOrigin == "x"}
            onIconClick={() => onSelectCurveOrigin("x")}
          />
          <SVGIcon
            sprite={yIcon}
            active={curveOrigin == "y"}
            onIconClick={() => onSelectCurveOrigin("y")}
          />
          <SVGIcon
            sprite={scaleIcon}
            active={scaleTo == "filtered"}
            onIconClick={() =>
              onSelectScaleTo(scaleTo == "total" ? "filtered" : "total")
            }
          />
        </MenuDisplaySimple>
        <MenuDisplaySimple name="maximum span">
          <div className={styles.full_height}>
            {meta && meta.assembly && meta.assembly.span < maxSpan && (
              <span
                className={styles.reset}
                onClick={() => onChangeMaxSpan(meta.assembly.span)}
              >
                reset
              </span>
            )}
            <NumericInput
              initialValue={maxSpan}
              minValue={meta.assembly ? meta.assembly.span : 0}
              onChange={onChangeMaxSpan}
            />
          </div>
        </MenuDisplaySimple>
        <MenuDisplaySimple name="maximum count">
          <div className={styles.full_height}>
            {records < maxCount && (
              <span
                className={styles.reset}
                onClick={() => onChangeMaxCount(records)}
              >
                reset
              </span>
            )}
            <NumericInput
              initialValue={maxCount}
              minValue={records ? records : 0}
              onChange={onChangeMaxCount}
            />
          </div>
        </MenuDisplaySimple>
        <ReferenceAssembliesList />
        {displayTotal}
      </span>
    );
    let snail = (
      <span>
        <MenuDisplaySimple name="snail origin">
          <SVGIcon
            sprite={invertIcon}
            active={snailOrigin == "center"}
            onIconClick={() =>
              onSelectSnailOrigin(snailOrigin == "center" ? "outer" : "center")
            }
          />
          <SVGIcon
            sprite={scaleIcon}
            active={scaleTo == "filtered"}
            onIconClick={() =>
              onSelectScaleTo(scaleTo == "total" ? "filtered" : "total")
            }
          />
        </MenuDisplaySimple>
      </span>
    );
    switch (view) {
      case "blob":
        context = blob;
        break;
      case "busco":
        break;
      case "cumulative":
        context = cumulative;
        break;
      case "report":
        context = (
          <span>
            {blob}
            {cumulative}
            {snail}
          </span>
        );
        break;
      case "snail":
        context = snail;
        break;
      case "table":
        break;
      case "treemap":
        break;
    }
    let showBlob;
    if (data.axes && data.axes.y.values.length > 0) {
      showBlob = true;
    }
    return (
      <div className={styles.menu_outer}>
        <MenuSwitchView />
        <MenuDataset
          key={datasetId}
          id={datasetId}
          active={false}
          onDatasetClick={() => {}}
          onDatasetMount={() => {}}
        />
        <MenuDisplaySimple invert={false}>
          <TextIcon
            title="interactive"
            active={!isStatic}
            onIconClick={() => onToggleStatic(view, datasetId, isStatic)}
          />
          {hasStatic && (
            <TextIcon
              title="static"
              active={isStatic}
              onIconClick={() => onToggleStatic(view, datasetId, isStatic)}
            />
          )}
        </MenuDisplaySimple>
        <MenuDisplaySimple invert={false}>
          {showBlob && (
            <TextIcon
              title="blob"
              active={view == "blob"}
              onIconClick={() => onSelectView("blob")}
            />
          )}
          <TextIcon
            title="busco"
            active={view == "busco"}
            onIconClick={() => onSelectView("busco")}
          />
          <TextIcon
            title="cumulative"
            active={view == "cumulative"}
            onIconClick={() => onSelectView("cumulative")}
          />
          <TextIcon
            title="detail"
            active={view == "detail"}
            onIconClick={() => onSelectView("detail")}
          />
          {isStatic || (
            <TextIcon
              title="report"
              active={view == "report"}
              onIconClick={() => onSelectView("report")}
            />
          )}
          <TextIcon
            title="snail"
            active={view == "snail"}
            onIconClick={() => onSelectView("snail")}
          />
          {isStatic || (
            <TextIcon
              title="table"
              active={view == "table"}
              onIconClick={() => onSelectView("table")}
            />
          )}
        </MenuDisplaySimple>
        <MenuDisplaySet name={view}>
          {(isStatic && view != "blob") || context}
        </MenuDisplaySet>
        {isStatic || <Palettes />}
        {isStatic || (
          <MenuDisplaySimple name="png resolution (px)">
            <NumericInput
              initialValue={pngResolution}
              onChange={onChangePngResolution}
            />
          </MenuDisplaySimple>
        )}
        <MenuDisplaySimple name="static threshold">
          <NumericInput
            initialValue={staticThreshold}
            onChange={(value) => {
              onChangeStaticThreshold(value);
              setTimeout(() => {
                if (value >= records && isStatic) {
                  onToggleStatic(view, datasetId, isStatic);
                } else if (value < records && !isStatic) {
                  onToggleStatic(view, datasetId, isStatic);
                }
              }, 500);
            }}
          />
        </MenuDisplaySimple>
        <MenuDisplaySimple name="nohit threshold">
          <NumericInput
            initialValue={nohitThreshold}
            onChange={(value) => {
              onChangeNohitThreshold(value);
            }}
          />
        </MenuDisplaySimple>
        <MenuDisplaySimple name="use larger fonts">
          <SVGIcon
            sprite={fontSizeIcon}
            active={largeFonts}
            onIconClick={() => onToggleLargeFonts(!largeFonts)}
          />
        </MenuDisplaySimple>

        <ToolTips set="settingsMenu" />
      </div>
    );
  }
}
// <br/>
// <label htmlFor='transform_x'>x position: {transform.x} </label>
// <br/>
// <input id='transform_x' onChange={(e)=>onChangeTransform({x:e.target.value})} type="range" value={transform.x} min="0" max="1000" step="50"/>
// <br/>
// <label htmlFor='transform_order'>order: {transform.order} </label>
// <br/>
// <input id='transform_order' onChange={(e)=>onChangeTransform({order:e.target.value})} type="range" value={transform.order} min="0.25" max="3" step="0.25"/>
// <br/>
// <label htmlFor='transform_factor'>factor: {transform.factor} </label>
// <br/>
// <input id='transform_factor' onChange={(e)=>onChangeTransform({factor:e.target.value})} type="range" value={transform.factor} min="-1" max="1" step="0.1"/>

class MenuDisplayMain extends React.Component {
  constructor(props) {
    super(props);
    this.resChange = null;
    this.mapDispatchToProps = (dispatch) => {
      return {
        onSelectShape: (shape) => dispatch(choosePlotShape(shape)),
        onChangeResolution: (value) => {
          this.props.clearTimeouts();
          this.props.setTimeout(
            () => dispatch(choosePlotResolution(value)),
            1000
          );
          return dispatch(setPlotResolution(value));
        },
        onChangeGraphics: (graphics) => dispatch(choosePlotGraphics(graphics)),
        onChangeThreshold: (threshold) =>
          dispatch(chooseSVGThreshold(threshold)),
        onChangePngResolution: (resolution) =>
          dispatch(choosePngResolution(resolution)),
        onChangeStaticThreshold: (threshold) =>
          dispatch(chooseStaticThreshold(threshold)),
        onChangeNohitThreshold: (threshold) =>
          dispatch(chooseNohitThreshold(threshold)),
        onChangeCircleLimit: (limit) => dispatch(chooseCircleLimit(limit)),
        onSelectReducer: (reducer) => dispatch(chooseZReducer(reducer)),
        onSelectScale: (scale) => dispatch(chooseZScale(scale)),
        onChangePlotScale: (plotScale) => dispatch(choosePlotScale(plotScale)),
        onSelectView: (view) => dispatch(chooseView(view)),
        onToggleStatic: (view, datasetId, isStatic) =>
          dispatch(toggleStatic(view, datasetId, isStatic)),
        // onToggleAdjustCoverage: (bool) => dispatch(chooseAdjustCoverage(bool)),
        onToggleLargeFonts: (bool) => dispatch(chooseLargeFonts(bool)),
        onChangeErrorBars: (value) => dispatch(chooseErrorBars(value)),
        onChangeWindowSize: (value) => dispatch(chooseWindowSize(value)),
        onSelectCurveOrigin: (origin) => dispatch(chooseCurveOrigin(origin)),
        onSelectScaleTo: (origin) => dispatch(chooseScaleTo(origin)),
        onSelectSnailOrigin: (origin) => dispatch(chooseSnailOrigin(origin)),
        onChangeTransform: (object) => dispatch(setTransformFunction(object)),
        onChangeShowTotal: (bool) => dispatch(chooseShowTotal(bool)),
        onChangeSideMax: (value) => dispatch(chooseSideMax(value)),
        onChangeMaxSpan: (value, minValue) => {
          if (minValue) {
            value = Math.max(value, minValue);
          }
          dispatch(chooseMaxSpan(value));
        },
        onChangeMaxCount: (value, minValue) => {
          if (minValue) {
            value = Math.max(value, minValue);
          }
          dispatch(chooseMaxCount(value));
        },
        onChangeAxisRange: (values, remove) =>
          dispatch(queryToStore({ values, remove, action: "FILTER" })),
        changeQueryParams: (obj) =>
          dispatch(
            setQueryString(
              Object.keys(obj)
                .map((k) => `${k}=${obj[k]}`)
                .join("&")
            )
          ),
      };
    };

    this.mapStateToProps = (state, props) => {
      return {
        title: getAxisTitle(state, "z"),
        records: getRecordCount(state),
        meta: getSelectedDatasetMeta(state),
        fields: getFields(state),
        shape: getPlotShape(state),
        resolution: getPlotResolution(state),
        graphics: getPlotGraphics(state),
        threshold: getSVGThreshold(state),
        pngResolution: getPngResolution(state),
        staticThreshold: getStaticThreshold(state),
        nohitThreshold: getNohitThreshold(state),
        circleLimit: getCircleLimit(state),
        reducer: getZReducer(state),
        scale: getZScale(state),
        plotScale: getPlotScale(state),
        curveOrigin: getCurveOrigin(state),
        scaleTo: getScaleTo(state),
        snailOrigin: getSnailOrigin(state),
        transform: getTransformFunctionParams(state),
        view: getView(state),
        isStatic: getStatic(state),
        datasetId: getDatasetID(state),
        hasStatic: getStaticFields(state),
        busco: getBuscoSets(state),
        data: getMainPlotData(state),
        parsed: getParsedQueryString(state),
        showTotal: getShowTotal(state),
        sideMax: getSideMax(state),
        maxSpan: getMaxSpan(state),
        maxCount: getMaxCount(state),
        // adjustCoverage: getAdjustCoverage(state),
        largeFonts: getLargeFonts(state),
        errorBars: getErrorBars(state),
        windowSize: getWindowSize(state),
      };
    };
  }

  render() {
    const DisplayMain = connect(
      this.mapStateToProps,
      this.mapDispatchToProps
    )(DisplayMenu);
    return <DisplayMain {...this.props} />;
  }
}
export default Timeout(MenuDisplayMain);
