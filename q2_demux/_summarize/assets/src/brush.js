import * as d3 from 'd3';

export default function addBrush(svg, plot, props, x, y, x0, y0, xAxis, yAxis) {
  const brush = d3.brush();
  let idleTimeout;
  const idleDelay = 350;

  const idled = () => {
    idleTimeout = null;
  };

  const zoom = () => {
    const t = svg.transition().duration(750);
    svg.select('.axis--x').transition(t).call(xAxis);
    svg.select('.axis--y').transition(t).call(yAxis);
    plot
      .domain(y.domain())
      .width((x.range()[1] - x.range()[0]) / (x.domain()[1] - x.domain()[0]) / 2);
    svg.selectAll('.boxplot').transition(t)
      .call(plot)
      .attr('transform', d => `translate(${x(d[0])}, ${props.margin.top})`);
  };

  const brushEnded = () => {
    const s = d3.event.selection;
    if (!s) {
      if (!idleTimeout) {
        idleTimeout = setTimeout(idled, idleDelay);
        return idleTimeout;
      }
      x.domain(x0);
      y.domain(y0);
    } else {
      x.domain([s[0][0], s[1][0]].map(x.invert, x));
      y.domain([s[1][1], s[0][1]].map(y.invert, y));
      svg.select('.brush').call(brush.move, null);
    }
    return zoom();
  };

  brush.on('end', brushEnded);

  svg.append('g')
      .attr('class', 'brush')
      .call(brush);
}
