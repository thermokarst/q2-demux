// Heavily modified variant of http://bl.ocks.org/jensgrubert/7789216
// which is based off of https://gist.github.com/mbostock/4061502
// which is written and marked as being under GPL-3

export default function addBox(d3) {
  /* eslint-disable no-shadow */
  // Inspired by http://informationandvisualization.de/blog/box-plot
  function boxWhiskers(d) {
    return [0, d.length - 1];
  }

  function boxQuartiles(d) {
    return [
      d3.quantile(d, 0.25),
      d3.quantile(d, 0.5),
      d3.quantile(d, 0.75),
    ];
  }

  d3.box = () => {  // eslint-disable-line no-param-reassign
    let width = 1;
    let height = 1;
    let duration = 750;
    let domain = null;
    let value = Number;
    let whiskers = boxWhiskers;
    let quartiles = boxQuartiles;
    let tickFormat = null;

    function box(gs) {
      gs.each(function draw(data, i) {
        const d = data[1].sort(d3.ascending);

        const g = d3.select(this);
        const n = d.length;
        const min = d[0];
        const max = d[n - 1];

        d.quartiles = quartiles(d);
        const quartileData = d.quartiles;

        const whiskerIndices = whiskers && whiskers.call(this, d, i);
        const whiskerData = whiskerIndices && whiskerIndices.map(wi => d[wi]);

        const outlierIndices = whiskerIndices
          ? d3.range(0, whiskerIndices[0]).concat(d3.range(whiskerIndices[1] + 1, n))
          : d3.range(n);

        const x1 = d3.scaleLinear()
          .domain((domain && domain.call(this, d, i)) || [min, max])
          .range([height, 0]);

        const x0 = this.existingChart || d3.scaleLinear()
          .domain([0, Infinity])
          .range(x1.range());

        this.existingChart = x1;

        const center = g.selectAll('line.center')
          .data(whiskerData ? [whiskerData] : []);
        center.enter().append('line', 'rect')
          .attr('class', 'center')
          .attr('x1', 0)
          .attr('y1', d => x0(d[0]))
          .attr('x2', 0)
          .attr('y2', d => x0(d[1]))
          .style('opacity', 1e-6)
          .attr('stroke-dasharray', '3,3')
          .attr('stroke', 'black')
          .attr('stroke-width', '1px')
        .transition()
          .duration(duration)
          .style('opacity', 1)
          .attr('y1', d => x1(d[0]))
          .attr('y2', d => x1(d[1]));

        center.transition()
          .duration(duration)
          .style('opacity', 1)
          .attr('x1', 0)
          .attr('x2', 0)
          .attr('y1', d => x1(d[0]))
          .attr('y2', d => x1(d[1]));

        center.exit().transition()
          .duration(duration)
          .style('opacity', 1e-6)
          .attr('y1', d => x1(d[0]))
          .attr('y2', d => x1(d[1]))
          .remove();

      // Update innerquartile box.
        const iqBox = g.selectAll('rect.box')
          .data([quartileData]);

        iqBox.enter().append('rect')
          .attr('class', 'box')
          .attr('x', -width / 2)
          .attr('y', d => x0(d[2]))
          .attr('width', width)
          .attr('height', d => x0(d[0]) - x0(d[2]))
          .attr('fill', 'steelblue')
          .attr('stroke', 'black')
          .attr('stroke-width', '1px')
          .on('mouseover', function mouseover() {
            d3.select(this)
              .attr('fill', 'skyblue');
          })
          .on('mouseout', function mouseout() {
            d3.select(this)
              .attr('fill', 'steelblue');
          })
        .transition()
          .duration(duration)
          .attr('y', d => x1(d[2]))
          .attr('height', d => x1(d[0]) - x1(d[2]));

        iqBox.transition()
          .duration(duration)
          .attr('x', -width / 2)
          .attr('y', d => x1(d[2]))
          .attr('width', width)
          .attr('height', d => x1(d[0]) - x1(d[2]));

      // Update median line.
        const medianLine = g.selectAll('line.median')
          .data([quartileData[1]]);

        medianLine.enter().append('line')
          .attr('class', 'median')
          .attr('x1', -width / 2)
          .attr('y1', x0)
          .attr('x2', width / 2)
          .attr('y2', x0)
          .attr('stroke', 'black')
          .attr('stroke-width', '1px')
        .transition()
          .duration(duration)
          .attr('y1', x1)
          .attr('y2', x1);

        medianLine.transition()
          .duration(duration)
          .attr('x1', -width / 2)
          .attr('x2', width / 2)
          .attr('y1', x1)
          .attr('y2', x1);

      // Update whiskers.
        const whisker = g.selectAll('line.whisker')
          .data(whiskerData || []);

        whisker.enter().insert('line', 'circle, text')
          .attr('class', 'whisker')
          .attr('x1', -width / 2)
          .attr('y1', x0)
          .attr('x2', width / 2)
          .attr('y2', x0)
          .style('opacity', 1e-6)
          .attr('stroke', 'black')
          .attr('stroke-width', '1px')
        .transition()
          .duration(duration)
          .attr('y1', x1)
          .attr('y2', x1)
          .style('opacity', 1);

        whisker.transition()
          .duration(duration)
          .attr('x1', -width / 2)
          .attr('x2', width / 2)
          .attr('y1', x1)
          .attr('y2', x1)
          .style('opacity', 1);

        whisker.exit().transition()
          .duration(duration)
          .attr('y1', x1)
          .attr('y2', x1)
          .style('opacity', 1e-6)
          .remove();

        const outlier = g.selectAll('circle.outlier')
          .data(outlierIndices, Number);

        outlier.enter().insert('circle', 'text')
          .attr('class', 'outlier')
          .attr('r', width / 4)
          .attr('cx', 0)
          .attr('cy', i => x0(d[i]))
          .style('opacity', 1e-6)
          .attr('fill', 'none')
          .attr('stroke', 'black')
        .transition()
          .duration(duration)
          .attr('cy', i => x1(d[i]))
          .style('opacity', 1);

        outlier.transition()
          .duration(duration)
          .attr('cy', i => x1(d[i]))
          .attr('r', width / 8)
          .attr('cx', 0)
          .style('opacity', 1);

        outlier.exit().transition()
          .duration(duration)
          .attr('cy', i => x1(d[i]))
          .style('opacity', 1e-6)
          .remove();


        const format = tickFormat || x1.tickFormat(8);

        const boxTick = g.selectAll('text.box')
          .data(quartileData);

        boxTick.transition()
          .duration(duration)
          .text(format)
          .attr('y', x1);

        const whiskerTick = g.selectAll('text.whisker')
          .data(whiskerData || []);

        whiskerTick.transition()
          .duration(duration)
          .text(format)
          .attr('y', x1)
          .style('opacity', 1);

        whiskerTick.exit().transition()
          .duration(duration)
          .attr('y', x1)
          .style('opacity', 1e-6)
          .remove();

        d3.select(this).on('mouseover', function mouseover(d) {
          const svg = d3.select(this.parentNode).node();
          const plotContainer = d3.select(svg.parentNode);
          plotContainer
          .select('.stats')
          .select('tbody')
          .selectAll('tr')
            .data([
              ['Position Number', d[0]],
              ['Minimum', whiskerData[0]],
              ['1st Quartile', quartileData[0]],
              ['Median', quartileData[1]],
              ['3rd Quartile', quartileData[2]],
              ['Maximum', whiskerData[1]],
            ])
          .selectAll('td')
            .data(d => d)
            .text(d => d);
        });
      });
    }

    box.width = (x) => {
      if (!arguments.length) return width;
      width = x;
      return box;
    };

    box.height = (x) => {
      if (!arguments.length) return height;
      height = x;
      return box;
    };

    box.tickFormat = (x) => {
      if (!arguments.length) return tickFormat;
      tickFormat = x;
      return box;
    };

    box.duration = (x) => {
      if (!arguments.length) return duration;
      duration = x;
      return box;
    };

    const functor = x => () => x;

    box.domain = (x) => {
      if (!arguments.length) return domain;
      domain = x == null ? x : functor(x);
      return box;
    };

    box.value = (x) => {
      if (!arguments.length) return value;
      value = x;
      return box;
    };

    box.whiskers = (x) => {
      if (!arguments.length) return whiskers;
      whiskers = x;
      return box;
    };

    box.quartiles = (x) => {
      if (!arguments.length) return quartiles;
      quartiles = x;
      return box;
    };

    return box;
  };
  /* eslint-enable no-shadow */
}
