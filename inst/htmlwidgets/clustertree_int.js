HTMLWidgets.widget({

  name: 'clustertree_int',

  type: 'output',

  factory: function(el, width, height) {

    // TODO: define shared variables for this instance

    return {

      renderValue: function(x) {

        // TODO: code to render the widget, e.g.
        //el.innerText = x.message;
        var scatter_div = document.createElement("div");
        scatter_div.setAttribute("id", "scatter_div");
        scatter_div.innerHTML = x.svg_other;

        var ctree_div = document.createElement("div");
        ctree_div.setAttribute("id", "ctree_div");
        ctree_div.innerHTML = x.svg_ctree;

        var input_rng = document.createElement("input");
        input_rng.setAttribute("type", "range");
        input_rng.setAttribute("min", x.tree_rng.min);
        input_rng.setAttribute("max", x.tree_rng.max);
        input_rng.setAttribute("value", x.tree_rng.min);
        input_rng.setAttribute("value", x.tree_rng.step);

        // Add both to plot
        el.appendChild(ctree_div);
        el.appendChild(scatter_div);
        el.appendChild(input_rng);


        // use this to sort of make our diagram responsive
        //  or at a minimum fit within the bounds set by htmlwidgets
        //  for the parent container
        function makeResponsive(el){
          var svg1 = document.getElementById("scatter_div").getElementsByTagName('svg')[0];
          var svg2 = document.getElementById("ctree_div").getElementsByTagName('svg')[0];
          //var document.getElementById("scatter_div");
          if(svg1){
            if(svg1.width) {svg1.removeAttribute("width")};
            if(svg1.height) {svg1.removeAttribute("height")};
            svg1.style.width = "50%";
            svg1.style.height = "50%";
          }
          if(svg2){
            if(svg2.width) {svg2.removeAttribute("width")};
            if(svg2.height) {svg2.removeAttribute("height")};
            svg2.style.width = "50%";
            svg2.style.height = "50%";
          }
        };

        makeResponsive(el);
      },

      setCutLine : function(params) {
        var line = document.getElementById("cut_line");

      },

      resize: function(width, height) {

        // TODO: code to re-render the widget with a new size

      }

    };
  }
});