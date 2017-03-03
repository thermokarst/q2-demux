// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';
import box from './box';


const plot = (data, svg, props) => {
  const maxX = d3.max(data, d => d[0]);
  const x0 = [0, maxX];
  const y0 = [0, 100];
  const x = d3.scaleLinear().domain(x0).range([props.margin.left, props.width]);
  // const x = d3.scaleBand().domain(data.map(d => d[0])).rangeRound([0, props.width])
  const y = d3.scaleLinear().domain(y0).range([props.height - props.margin.bottom, props.margin.top]);

  var xAxis = d3.axisTop(x).ticks(12),
      yAxis = d3.axisRight(y).ticks(12 * props.height / props.width);

  const brush = d3.brush().on("end", brushEnded);
  let idleTimeout;
  const idleDelay = 350;

  var chart = d3.box()
    .whiskers(iqr(1.5))
    .height(props.height)
    .domain([0, 100])
    .width(0.7)
    .showLabels(false);

  function iqr(k) {
    return function(d, i) {
      var q1 = d.quartiles[0],
          q3 = d.quartiles[2],
          iqr = (q3 - q1) * k,
          i = -1,
          j = d.length;
      while (d[++i] < q1 - iqr);
      while (d[--j] > q3 + iqr);
      return [i, j];
    }
  }

  // svg.selectAll("circle")
  //   .data(data)
  //   .enter().append("circle")
  //     .attr("cx", d => x(d.x))
  //     .attr("cy", d => y(d.y))
  //     .attr("r", 2.5)
  //     .attr("fill", d => d3.interpolateViridis(d.y / 100))
  //     .attr('fill-opacity', 0.5)
  //     .call(chart);
  svg.selectAll(".box")
    .data(data)
  .enter().append("g")
    .attr("transform", d => `translate(${x(d[0])}, ${props.margin.top})`)
    .call(chart);

  svg.append("g")
      .attr("class", "axis axis--x")
      .attr("transform", "translate(0," + (props.height - 10) + ")")
      .call(xAxis);

  svg.append("g")
      .attr("class", "axis axis--y")
      .attr("transform", "translate(10,0)")
      .call(yAxis);

  svg.selectAll(".domain")
      .style("display", "none");

  svg.append("g")
      .attr("class", "brush")
      .call(brush);

  function brushEnded() {
    var s = d3.event.selection;
    if (!s) {
      if (!idleTimeout) return idleTimeout = setTimeout(idled, idleDelay);
      x.domain(x0);
      y.domain(y0);
    } else {
      x.domain([s[0][0], s[1][0]].map(x.invert, x));
      y.domain([s[1][1], s[0][1]].map(y.invert, y));
      svg.select(".brush").call(brush.move, null);
    }
    zoom();
  }

  const idled = () => {
    idleTimeout = null;
  }

  const zoom = () => {
    var t = svg.transition().duration(750);
    svg.select(".axis--x").transition(t).call(xAxis);
    svg.select(".axis--y").transition(t).call(yAxis);
    svg.selectAll(".box").transition(t)
      .attr("transform", d => `translate(${x(d[0])}, ${props.margin.top})`)
  }
}

const initializePlot = (data) => {
  box(d3);
  const margin = { top: 10, right: 30, bottom: 30, left: 30 };
  const width = d3.select('#chartContainer').node().offsetWidth;
  const props = {
    margin,
    width: width - margin.left - margin.right,
    height: ((width * 9) / 16) - margin.top - margin.bottom };

  const svg = d3
    .select('#chartContainer')
  .append('svg')
    .attr('width', props.width + props.margin.left + props.margin.right)
    .attr('height', props.height + props.margin.top + props.margin.bottom);

  const maxLen = d3.max(data, d => d.length);
  const processedData = new Array(maxLen);
  for (const point of data) {
    for (let i = 0; i < point.length; i++) {
      processedData[i] = processedData[i] || new Array(2);
      if (!processedData[i][0]) {
        processedData[i][0] = i
      }
      if (!processedData[i][1]){
        processedData[i][1] = []
      }
      processedData[i][1].push(point[i]);
    }
  }

  plot(processedData, svg, props);
}

export default initializePlot;
