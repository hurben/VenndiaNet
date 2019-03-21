// by @lazykuna
// MIT license
//
// refered Chord Diagram: https://bl.ocks.org/mbostock/4062006
//

/*
option
innerRadius
outerRadius
padAngle
*/

d3.multichord = function(option) {
  var _outerPadding = (option&&option.outerPadding?option.outerPadding:40),    // outer padding
      _innerPadding = (option&&option.innerPadding?option.innerPadding:80),    // inner padding
      _padAngle = ((option&&option.padAngle)?option.padAngle:0);

  var multichord = function(svg, graph) {
    // check svg element
    //
    if (svg.node().tagName !== 'svg')
    {
      let parentWidth = svg.select('svg').node().parentNode.clientWidth;
      let parentHeight = svg.select('svg').node().parentNode.clientHeight;
      svg = svg.select('svg')
      svg
        .attr('width', parentWidth)
        .attr('height', parentHeight);
    }

    // initialize
    var width = +svg.attr("width"),
        height = +svg.attr("height"),
        outerRadius = Math.min(width, height) * 0.5 - _outerPadding,   // outer radius
        innerRadius = outerRadius - _innerPadding;     // thickness of inner round

    console.assert(outerRadius > 0);
    console.assert(innerRadius > 0);

    var padAngle = _padAngle;
    var formatValue = d3.formatPrefix(",.0", 1e3);

    var arc = d3.arc()
        .innerRadius(innerRadius)
        .outerRadius(outerRadius);

    /*
    var color = d3.scaleOrdinal()
        .domain(d3.range(4))
        .range(["#000000", "#FFDD89", "#957244", "#F26223"]);
    */
    var color = d3.scaleOrdinal(d3.schemeCategory20); //d3.scale.category20();


    var cos = Math.cos;
    var sin = Math.sin;
    var pi = Math.PI;
    var tau = pi * 2;
    var max = Math.max;
    var halfPi = pi / 2;
    var tooltip = d3.select("body")
      .append("div")
      .attr("class", "tooltip-multichord")
      .style("z-index", "10")
      .style("visibility", "hidden")
      .style("position", "absolute");


    // event handler
    var onGroupMouseOver = function(d,i)
    {
      // collect all related sets
      var groups = [d.name];
      d3.select(".rbn-"+d.name).each(function(d2,i){
        if (d2.name.indexOf("__") >= 0) {
          groups.push(d2.name);
        }
      });
      var groups_str = groups.join(",");
      // transition on all related sets
      d3.select(this).transition()
        .ease(d3.easeQuad)
        .duration("200")
        .style("opacity", 1);
      d3.selectAll(".rbn-"+d.name).transition()
        .ease(d3.easeQuad)
        .duration("200")
        .style("opacity", 1);
      // tooltip
      tooltip.html("(Union Set)<br>"+d.name);
      tooltip.style("visibility", "visible")
        .style("left", (d3.event.pageX + 10) + "px")   
        .style("top", (d3.event.pageY - 28) + "px");  
      // customized event (object, index, selected sets)
      if (option&&option.onGroupMouseOver)
      {
        option.onGroupMouseOver(d,i,groups);
      }
    };
    var onGroupMouseOut = function(d,i)
    {
      d3.select(this).transition()
        .ease(d3.easeQuad)
        .duration("200")
        .style("opacity", 0.5);
      d3.selectAll(".rbn-"+d.name).transition()
        .ease(d3.easeQuad)
        .duration("200")
        .style("opacity", 0.5);
      tooltip.style("visibility", "hidden");
      // customized event (object, index, selected sets)
      if (option&&option.onGroupMouseOut)
      {
        option.onGroupMouseOut(d,i);
      }
    };
    var onGroupClick = function(d,i)
    {
      var groups = [d.name];
      d3.selectAll(".rbn-"+d.name).each(function(d2,i){
        if (d2.name.indexOf("__") >= 0) {
          groups.push(d2.name);
        }
      });
      if (option&&option.onGroupClick)
      {
        option.onGroupClick(d,i,groups);
      }
    };
    var onRibbonMouseOver = function(d,i)
    {
      d3.select(this).transition()
        .ease(d3.easeQuad)
        .duration("200")
        .style("opacity", 1);
      var names = d.name.split('__');
      tooltip.html(""+names.join('<br> & '));
      tooltip.style("visibility", "visible")
        .style("left", (d3.event.pageX + 10) + "px")   
        .style("top", (d3.event.pageY - 28) + "px");  
      // customized event (object, index, selected sets)
      if (option&&option.onGroupMouseOver)
      {
        option.onGroupMouseOver(d,i,[d.name]);
      }
    };
    var onRibbonMouseOut = function(d,i)
    {
      d3.select(this).transition()
        .ease(d3.easeQuad)
        .duration("200")
        .style("opacity", 0.5);
      tooltip.style("visibility", "hidden");
      // customized event (object, index, selected sets)
      if (option&&option.onGroupMouseOut)
      {
        option.onGroupMouseOut(d,i);
      }
    };
    var onRibbonClick = function(d,i)
    {
      if (option&&option.onGroupClick)
      {
        option.onGroupClick(d,i,[d.name]);
      }
    };


    // converts data format: graph --> sets(chords)
    function convertToSet(graph)
    {
      var sets = {};
      if (!graph.groups)
      {
        // make it into dict
        for (var i=0; i<graph.nodes.length; i++)
        {
          var node = graph.nodes[i];
          if (!(node.group in sets))
          {
            sets[node.group] = [];
            sets[node.group].color = color(node.group);
            sets[node.group].count = 0;
          }
          // make color scheme
          if (!('color' in node)) {
            //node.color = color(node.group);
          }
          sets[node.group].push(node);
          sets[node.group].count++;
        }
      } else {
        for (var i=0; i<graph.groups.length; i++)
        {
          sets[graph.groups[i].label] = graph.groups[i];
        }
      }
      return sets;
    }

    function getRibbonPath(ribbons)
    {
      var context = d3.path();

      for (var i=0; i<ribbons.length; i++)
      {
        var r = innerRadius,
          sa = ribbons[i][0] - halfPi,
          ea = ribbons[i][1] - halfPi,
          sx = r * cos(sa),
          sy = r * sin(sa);
        if (i==0)
        {
          // starting point using last set's angle
          var ea2 = ribbons[ribbons.length-1][1] - halfPi;
          var sx1 = r * cos(ea2);
          var sy1 = r * sin(ea2);
          context.moveTo(sx1, sy1);
          // draw outer curve 
          // then starting point is automatically set
          context.quadraticCurveTo(0, 0, sx, sy);
        }
        else
        {
          // draw curve between set
          context.quadraticCurveTo(0, 0, sx, sy);
        }
        // draw arc
        // TODO: (if and only if angle is different)
        context.arc(0, 0, r, sa, ea);
      }
      context.closePath();

      return context.toString();
    }


    // set, angle calculator
    function multiset_chord(sets)
    {
      // parameter = padding_angle
      // 
      var padding_angle = padAngle;

      // 1. calculate angle for sets
      // 1-1. count total element count first
      //
      var chorddict = {};   // key: [startangle, endangle, value]
      var chordcount_total = 0;
      var ccount_total = 0;
      var idx = 0;
      for (var n in sets)
      {
        var chord_intersection = n.split('__');
        var chord_count = sets[n].count;
        for (var j=0; j<chord_intersection.length; j++)
        {
          var name = chord_intersection[j];
          if(!(name  in chorddict)) {
            chorddict[name] = {'startAngle': 0, 'endAngle':0, 'value':0, 'curAngle':0, 'index': idx, 'subindex': 0, 'color': sets[n].color};
            idx++;
          }
          chorddict[name]['value'] += chord_count;
          chordcount_total += chord_count;
        }
        ccount_total += chord_count;
      }

      // 1-2. calculate angle using count
      //
      var chordangle_start = 0;
      var chordangle_total = Math.PI*2 - padding_angle*(Object.keys(chorddict).length);
      for (var n in chorddict)
      {
        chorddict[n]['startAngle']
          = chorddict[n]['curAngle']
          = chordangle_start;
        var chordangle = chordangle_total * chorddict[n]['value'] / chordcount_total;
        chorddict[n]['endAngle'] = chordangle+chordangle_start;
        chordangle_start += chordangle + padding_angle;
      }
      
      // 2. calculate angle for ribbons
      //
      var ribbons = [];
      for (var n in sets)
      {
        var chord_intersection = n.split('__');
        var chord_count = sets[n].count;
        var delta_angle = chordangle_total * chord_count / chordcount_total;
        var ribbon = [];
        var clr = chorddict[ chord_intersection[0] ]['color'];   // color scheme of ribbon
        /*
        if (chord_intersection.length > 1) {
          clr += chord_intersection.length + 10; // TODO: (color scheme) need to be more solid?
        }*/
        for (var j=0; j<chord_intersection.length; j++)
        {
          var curname = chord_intersection[j];
          var sa = chorddict[curname]['curAngle'];
          chorddict[curname]['curAngle']+=delta_angle;
          var ea = chorddict[curname]['curAngle'];
          ribbon.push(
            [sa, ea]
          );
        }
        ribbon.value = 0;
        //ribbon.index = curindex;
        ribbon.color = sets[n].color;
        ribbon.sets = chord_intersection;
        ribbon.name = n;
        ribbons.push(ribbon);
      }

      // 3. summarize
      //
      var chordset = [];
      for (var n in chorddict)
      {
        chordset.push(Object.assign({}, chorddict[n], {'name': n}));
      }
      ribbons.groups = chordset;
      return ribbons;
    }

    // translate data
    var sets = convertToSet(graph);
    var chords = multiset_chord(sets);


    var g = svg.append("g")
        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
        .datum(chords);

    var arc = d3.arc()
        .innerRadius(innerRadius)
        .outerRadius(outerRadius);

    // render groups
    // 
    var group = g.append("g")
        .attr("class", "groups")
      .selectAll("g")
      .data(function(chords) { return chords.groups; })
      .enter().append("g");
    group.append("path")
        .style("fill", function(d) { if (d.color) return d.color; else return color(d.name); })
        .style("stroke", function(d) { if (d.color) return d3.rgb(d.color).darker(); else return d3.rgb(color(d.name)).darker(); })
        .style("opacity", 0.5)
        .attr("d", arc)
        .attr("class", "multichord-group")
    // onhover(groups)
      .on("mouseover", onGroupMouseOver)
      .on("mouseout", onGroupMouseOut)
      .on("click", onGroupClick);

    // render ribbons
    //
    g.append("g")
        .attr("class", "ribbons")
      .selectAll("path")
      .data(function(chords) { return chords; })
      .enter().append("path")
        .attr("d", getRibbonPath)
        .attr("class", function(d) { var cls="multichord-ribbon "; if (d.sets.length>1) { cls+="multiset "; } for (var s in d.sets) { cls += "rbn-"+d.sets[s]+" "; } return cls; })
        .style("opacity", 0.5)
        .style("fill", function(d) { if (d.color) return d.color; else return color(d.name); })
        .style("stroke", function(d) { if (d.color) return d3.rgb(d.color).darker(); else return d3.rgb(color(d.name)).darker(); })
    // onhover(ribbon)
      .on("mouseover", onRibbonMouseOver)
      .on("mouseout", onRibbonMouseOut)
      .on("click", onRibbonClick);


    return multichord;
  }

  multichord.innerPadding = function(v)
  {
    if (!v) { return _innerPadding; }
    _innerPadding = v;
    return multichord;
  }

  multichord.outerRadius = function(v)
  {
    if (!v) { return _outerPadding; }
    _outerPadding = v;
    return multichord;
  }

  multichord.padAngle = function(v)
  {
    if (!v) { return padAngle; }
    padAngle = v;
    return multichord;
  }

  multichord.onGroupMouseOver = function(v)
  {
    if (!v) { return option.onGroupMouseOver; }
    option.onGroupMouseOver = v;
    return multichord;
  }

  multichord.onGroupMouseOut = function(v)
  {
    if (!v) { return option.onGroupMouseOut; }
    option.onGroupMouseOut = v;
    return multichord;
  }

  multichord.onGroupClick = function(v)
  {
    if (!v) { return option.onGroupClick; }
    option.onGroupClick = v;
    return multichord;
  }

  return multichord;
};




// utility function
//
function createV4MultiChord(svg, graph, option)
{
  var multichord = d3.multichord(option);
  multichord(svg, graph);
  return multichord;
}




/*
var sets = {
  "A": {'length': 120},
  "B": {'length': 90},
  "C": {'length': 200},
  "A__C": {'length': 80},
}
*/
/*
var groupTick = group.selectAll(".group-tick")
  .data(function(d) { return groupTicks(d, 1e3); })
  .enter().append("g")
    .attr("class", "group-tick")
    .attr("transform", function(d) { return "rotate(" + (d.angle * 180 / Math.PI - 90) + ") translate(" + outerRadius + ",0)"; });

groupTick.append("line")
    .attr("x2", 6);

groupTick
  .filter(function(d) { return d.value % 5e3 === 0; })
  .append("text")
    .attr("x", 8)
    .attr("dy", ".35em")
    .attr("transform", function(d) { return d.angle > Math.PI ? "rotate(180) translate(-16)" : null; })
    .style("text-anchor", function(d) { return d.angle > Math.PI ? "end" : null; })
    .text(function(d) { return formatValue(d.value); });
*/



// input data type: nodes, group(no or string)
// TODO: make 'selection' work
// TODO: click event (with shift / alt)
