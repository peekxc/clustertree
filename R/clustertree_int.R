#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
clustertree_int <- function(width = NULL, height = NULL, elementId = NULL) {

  data(iris)
  cl_tree <- hclust(dist(iris[, 1:2]), method = "single")
  hgt <- cl_tree$height
  step_sz <- min(Filter(function(x) x != 0, abs(hgt[1:length(hgt) - 1] - hgt[2:length(hgt)])))

  ## Add line
  svg_ctree <- svglite:::inlineSVG({ plot(cl_tree); abline(h=0, col="red") });
  #svg_ctree <- gsub(pattern = "<line ([ [:alnum:]\"=]*)>(\\s+<[\\/]/line>\\s+<[\\/]/svg>)", "<line id = \"cut_line\" \\1 \\2", x = svg_ctree)
  svg_other <- svglite:::inlineSVG(plot(iris[, 1:2]))

  # forward options using x
  x = list(
    svg_ctree = svg_ctree,
    svg_other = svg_other,
    tree_rng = list(min=min(cl_tree$height), max=max(cl_tree$height), step = step_sz)
  )
  #

  # create widget
  htmlwidgets::createWidget(
    name = 'clustertree_int',
    x,
    width = width,
    height = height,
    package = 'clustertree',
    elementId = elementId
  )
}

#' Shiny bindings for clustertree_int
#'
#' Output and render functions for using clustertree_int within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a clustertree_int
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name clustertree_int-shiny
#'
#' @export
clustertree_intOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'clustertree_int', width, height, package = 'clustertree')
}

#' @rdname clustertree_int-shiny
#' @export
renderClustertree_int <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, clustertree_intOutput, env, quoted = TRUE)
}
