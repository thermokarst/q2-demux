// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';
import box from './box';


const plot = (data, props, container, paired) => {
  const plotContainer = d3.select(container);

  const svg = plotContainer
    .append('svg')
      .attr('class', paired ? 'col-xs-12' : 'col-xs-12 col-lg-10')
      .style('display', 'block')
      .style('margin', '0 auto')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', props.height + props.margin.top + props.margin.bottom);

  const panel = plotContainer
    .append('div')
      .attr('class', paired ? 'col-xs-12': 'col-xs-12 col-lg-2')
    .append('div')
      .attr('class', 'panel panel-default')
  panel.append('div')
    .attr('class', 'panel-heading')
    .html('Statistical Summary')
  panel.append('div')
    .attr('class', 'stats')
    .html(`<table class="table" style="margin-bottom: 0;"><thead><tr><th class="col-xs-9" style="font-size: 10px;">Hover over a boxplot to learn more...</th><th class="col-xs-3"></th></thead><tbody><tr><td>Base #</td><td>...</td></tr><tr><td>Minimum</td><td>...</td></tr><tr><td>1st Quartile</td><td>...</td></tr><tr><td>Median</td><td>...</td></tr><tr><td>3rd Quartile</td><td>...</td></tr><tr><td>Maximum</td><td>...</td></tr></tbody></table>`)

  const maxX = d3.max(data, d => d[0]);
  const x0 = [0, maxX];
  const y0 = [0, 100];
  const x = d3.scaleLinear().domain(x0).range([props.margin.left, props.width]);
  const y = d3.scaleLinear().domain(y0).range([props.height - props.margin.bottom, props.margin.top]);

  const xAxis = d3.axisBottom(x).ticks(12);
  const yAxis = d3.axisLeft(y).ticks(12 * props.height / props.width);

  const brush = d3.brush().on("end", brushEnded);
  let idleTimeout;
  const idleDelay = 350;

  var chart = d3.box()
    .whiskers(iqr(1.5))
    .height(props.height - props.margin.bottom - props.margin.top)
    .domain([0, 100])
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

  svg.append("g")
      .attr("class", "brush")
      .attr('z-index', -1)
      .call(brush);

  svg.selectAll(".box")
    .attr('font', '10px sans-serif')
    .data(data)
  .enter().append("g")
    .attr('class', 'boxplot')
    .attr("transform", d => `translate(${x(d[0])}, ${props.margin.top})`)
    .call(chart.width((x.range()[1] - x.range()[0]) / (x.domain()[1] - x.domain()[0]) / 2));

  svg.append('rect')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', props.margin.bottom * 4)
      .attr('y', props.height - props.margin.bottom)
      .attr('fill', 'white')

  svg.append('rect')
      .attr('width', props.margin.left * 4)
      .attr('x', -props.margin.left * 3)
      .attr('height', props.height + 200)
      .attr('fill', 'white')

  svg.append('rect')
    .attr('width', props.margin.right + 1000)
    .attr('height', props.height + 200)
    .attr('x', props.width + 10)
    .attr('fill', 'white')

  svg.append("g")
      .attr("class", "axis axis--x")
      .attr("transform", `translate(0, ${(props.height - props.margin.bottom)})`)
      .call(xAxis);

  svg.append("g")
      .attr("class", "axis axis--y")
      .attr("transform", `translate(${props.margin.left},0)`)
      .call(yAxis);

  svg.selectAll(".domain")
      .style("display", "none");

  svg.append('text')
      .attr("transform", "rotate(-90)")
      .attr("x", 0 - (props.height / 2))
      .attr('dy', '1em')
      .attr('font-size', '10px')
      .style("text-anchor", "middle")
      .text("Quality Score");

  svg.append('text')
      .attr("x", props.width / 2)
      .attr('y', props.height)
      .attr('dy', '1em')
      .attr('font-size', '10px')
      .style("text-anchor", "middle")
      .text("Sequence Base");

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
    svg.selectAll(".boxplot").transition(t)
      .call(chart.domain(y.domain())
                 .width((x.range()[1] - x.range()[0]) / (x.domain()[1] - x.domain()[0]) / 2))
      .attr("transform", d => `translate(${x(d[0])}, ${props.margin.top})`)
  }
}

const initializePlot = (data) => {
  box(d3);
  const margin = { top: 10, right: 30, bottom: 30, left: 30 };
  const width = d3.select('#forwardContainer').node().offsetWidth;
  const props = {
    margin,
    width: width - margin.left - margin.right,
    height: ((width * 9) / 16) - margin.top - margin.bottom };

  const paired = Object.keys(data[0]).length == 2 ? true : false;

  for (let direction of Object.keys(data[0])) {
    const maxLen = d3.max(data, d => d[direction].length);
    const processedData = new Array(maxLen);
    for (const point of data) {
      for (let i = 0; i < point[direction].length; i++) {
        processedData[i] = processedData[i] || new Array(2);
        if (!processedData[i][0]) {
          processedData[i][0] = i
        }
        if (!processedData[i][1]){
          processedData[i][1] = []
        }
        processedData[i][1].push(point.forward[i]);
      }
    }

    plot(processedData, props, `#${direction}Container`, paired);
  }
}

export default initializePlot;
