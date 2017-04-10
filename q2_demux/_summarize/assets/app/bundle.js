var app=webpackJsonpapp([0],[function(t,e,r){"use strict";function a(t){return t&&t.__esModule?t:{default:t}}Object.defineProperty(e,"__esModule",{value:!0}),e.init=void 0;var n=r(1),i=a(n);e.init=function(t){var e=arguments.length>1&&void 0!==arguments[1]?arguments[1]:void 0,r={forward:Object.keys(t).map(function(e){return[+e+1,t[e]]})};e&&(r.reverse=Object.keys(e).map(function(t){return[+t+1,e[t]]})),(0,i.default)(r)}},function(t,e,r){"use strict";function a(t){return t&&t.__esModule?t:{default:t}}function n(t){if(t&&t.__esModule)return t;var e={};if(null!=t)for(var r in t)Object.prototype.hasOwnProperty.call(t,r)&&(e[r]=t[r]);return e.default=t,e}Object.defineProperty(e,"__esModule",{value:!0});var i=r(2),l=n(i),o=r(3),s=a(o),u=r(4),c=function(t,e,r){var a=l.select(r),n=a.append("svg").attr("class","col-xs-12").style("display","block").style("margin","0 auto").attr("width",e.width+e.margin.left+e.margin.right).attr("height",e.height+e.margin.top+e.margin.bottom),i=a.append("div").attr("class","col-xs-12").append("div").attr("class","panel panel-default");i.append("div").attr("class","panel-heading").html("Statistical Summary");var o=i.append("div").attr("class","stats").append("table").attr("class","table").style("margin-bottom","0");o.append("thead").append("tr").selectAll("th").data([["Hover over a boxplot to learn more...",9],["",3]]).enter().append("th").text(function(t){return t[0]}).attr("class",function(t){return"col-xs-"+t[1]}).style("font-size","10px");var c=[["Position Number","..."],["Minimum","..."],["1st Quartile","..."],["Median","..."],["3rd Quartile","..."],["Maximum","..."]],d=o.append("tbody").selectAll("tr").data(c).enter().append("tr");d.selectAll("td").data(function(t){return t}).enter().append("td").text(function(t){return t});var f=(l.max(t,function(t){return t[0]}),t.map(function(t){return t[0]})),p=[0,45],m=l.scaleBand().domain(f).range([e.margin.left,e.width]),h=l.scaleLinear().domain(p).range([e.height-e.margin.bottom,e.margin.top]),g=l.axisBottom(m),v=l.axisLeft(h).ticks(9);n.attr("height",e.height+e.margin.bottom+e.margin.top),n.attr("width",m.range()[0]+m.range()[1]),(0,u.addBrush)(n,t,e,m,h,f,p,g,v),(0,s.default)(n,t,e,m,h),n.append("rect").attr("width",e.width+e.margin.left+e.margin.right).attr("height",4*e.margin.bottom).attr("y",e.height-e.margin.bottom).attr("fill","white"),n.append("rect").attr("width",4*e.margin.left).attr("x",3*-e.margin.left).attr("height",e.height+200).attr("fill","white"),n.append("rect").attr("width",e.margin.right+1e3).attr("height",e.height+200).attr("x",e.width+10).attr("fill","white"),n.append("g").attr("class","axis axis--x").attr("transform","translate(0, "+(e.height-e.margin.bottom)+")").call(g),n.append("g").attr("class","axis axis--y").attr("transform","translate("+e.margin.left+",0)").call(v),n.append("text").attr("transform","rotate(-90)").attr("x",0-e.height/2).attr("dy","1em").attr("font-size","10px").style("text-anchor","middle").text("Quality Score"),n.append("text").attr("x",e.width/2).attr("y",e.height).attr("dy","1em").attr("font-size","10px").style("text-anchor","middle").text("Sequence Base"),(0,u.updateXTicks)(n,m)},d=function(t){var e={top:10,right:30,bottom:30,left:30},r=l.select("#forwardContainer").node().offsetWidth,a={margin:e,width:r-e.left-e.right,height:9*r/16-e.top-e.bottom};Object.keys(t).forEach(function(e){c(t[e],a,"#"+e+"Container")})};e.default=d},,function(t,e,r){"use strict";function a(t){if(t&&t.__esModule)return t;var e={};if(null!=t)for(var r in t)Object.prototype.hasOwnProperty.call(t,r)&&(e[r]=t[r]);return e.default=t,e}function n(t,e,r,a,n){var i=a.bandwidth()<2?1:a.bandwidth()/2,o=i/2,s=t.selectAll(".container").data(e).attr("font","10px sans-serif");s.exit().remove();var u=s.enter().append("g").attr("class","container"),c=s.merge(u).attr("transform",function(t){return"translate("+(a(t[0])?a(t[0]):-100)+", "+r.margin.top+")"}).attr("display",function(t){return a(t[0])?"initial":"none"}).on("mouseover",function(){var t=l.select(this).data(),e=t[0][1],r=[e["25%"],e["50%"],e["75%"]],a=[e.min,e.max],n=l.select(this.parentNode).node(),i=l.select(n.parentNode);i.select(".stats").select("tbody").selectAll("tr").data([["Position Number",t[0][0]],["Minimum",a[0]],["1st Quartile",r[0]],["Median",r[1]],["3rd Quartile",r[2]],["Maximum",a[1]]]).selectAll("td").data(function(t){return t}).text(function(t){return t})}),d=c.selectAll("line.center").data(function(t){return[t]});d.exit().remove();var f=d.enter().append("line");d.merge(f).attr("class","center").attr("x1",i).attr("y1",function(t){return n(t[1].min)}).attr("x2",i).attr("y2",function(t){return n(t[1].max)}).attr("stroke-dasharray","2,2").attr("stroke-width",1).attr("stroke","black");var p=c.selectAll("rect.box").data(function(t){return[t]});p.exit().remove();var m=p.enter().append("rect");p.merge(m).attr("class","box").attr("x",o).attr("y",function(t){return n(t[1]["75%"])}).attr("width",i).attr("height",function(t){return n(t[1]["25%"])-n(t[1]["75%"])}).attr("fill","steelblue").attr("stroke-width",1).attr("stroke","black").on("mouseover",function(){l.select(this).attr("fill","skyblue")}).on("mouseout",function(){l.select(this).attr("fill","steelblue")});var h=c.selectAll("line.median").data(function(t){return[t]});h.exit().remove();var g=h.enter().append("line");h.merge(g).attr("class","median").attr("x1",o).attr("y1",function(t){return n(t[1]["50%"])}).attr("x2",i+o).attr("y2",function(t){return n(t[1]["50%"])}).attr("stroke-width",1).attr("stroke","black");var v=c.selectAll("line.lower-whisker").data(function(t){return[t]});v.exit().remove();var x=v.enter().append("line");v.merge(x).attr("class","lower-whisker").attr("x1",o).attr("y1",function(t){return n(t[1].min)}).attr("x2",i+o).attr("y2",function(t){return n(t[1].min)}).attr("stroke-width",1).attr("stroke","black");var y=c.selectAll("line.upper-whisker").data(function(t){return[t]});y.exit().remove();var b=y.enter().append("line");y.merge(b).attr("class","upper-whisker").attr("x1",o).attr("y1",function(t){return n(t[1].max)}).attr("x2",i+o).attr("y2",function(t){return n(t[1].max)}).attr("stroke-width",1).attr("stroke","black")}Object.defineProperty(e,"__esModule",{value:!0}),e.default=n;var i=r(2),l=a(i)},function(t,e,r){"use strict";function a(t){return t&&t.__esModule?t:{default:t}}function n(t){if(t&&t.__esModule)return t;var e={};if(null!=t)for(var r in t)Object.prototype.hasOwnProperty.call(t,r)&&(e[r]=t[r]);return e.default=t,e}function i(t,e){var r=t.selectAll(".axis--x .tick"),a=e.domain().length;a>48?r.style("display",function(t){return t%12?"none":"initial"}):a>24?r.style("display",function(t){return t%2?"none":"initial"}):r.style("display","initial")}function l(t,e,r,a,n,l,o,u,d){var f=s.brush(),p=void 0,m=350,h=(t.selectAll(".axis--x .tick"),function(){p=null}),g=function(){t.select(".axis--x").call(u),t.select(".axis--y").call(d),(0,c.default)(t,e,r,a,n),i(t,a)},v=function(){var e=s.event.selection;if(e){var r=s.scaleQuantize().domain(a.range()).range(a.domain()),i=r(e[0][0]),u=r(e[1][0]);a.domain(Array(u-i+1).fill().map(function(t,e){return i+e})),n.domain([e[1][1],e[0][1]].map(n.invert,n)),t.select(".brush").call(f.move,null)}else{if(!p)return p=setTimeout(h,m);a.domain(l),n.domain(o)}return g()};f.on("end",v),t.append("g").attr("class","brush").call(f)}Object.defineProperty(e,"__esModule",{value:!0}),e.updateXTicks=i,e.addBrush=l;var o=r(2),s=n(o),u=r(3),c=a(u)}]);