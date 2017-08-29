#' @title Plot a spanning tree
#' @name spanplot
#' @description Allows for visualizing a clustertree object as a spanning tree.
#' @param x a matrix to plot as a scatter plot.
#' @param which which hierarchy to plot. If RSL was used or k was sufficiently small, this may be only 1 hierarchy.
#' @param h an optional height to cut the tree at. Colors the corresponding connected components.
#' @param f optional function to plot a 2-dimensional representation of the connected components. See details.
#' @param ... arguments to pass to f. Otherwise unused.
#' If \code{f} is unspecified, then a simple two dimensional scatter plot is used to plot the nodes of the
#' spanning tree based on the following rules:
#' 1) if the original data is 1-dimensional, it's reflected symmetrically across the x and y axis (plotting on both dimensions)
#' 2) if the original data is 2-dimensional, it's plotted as a simple scatter plot (first column is X dimension).
#' 3) if the original data is n-dimensional, where n > 2, then the original data is mapped to the first two principle components.
#' If \code{f} is specified, it must take as its first argument a set of indices corresponding to the points in the
#' original data set, and it must return a set of two-dimensional points for those indices. These coordinates need not
#' be in the original data space. See the examples section for a demonstration of this.
#' @importFrom methods is
#' @import stats
#' @import graphics
#' @importFrom grDevices rainbow
#' @export
spanplot <- function(x, which = 1, h = NULL, f= NULL, ... ){
  if (!any(c("clustertree") %in% class(x))) stop("spanplot expects a 'clustertree' object")
  if (is.null(x[["mst"]])) stop("Cannot plot spanning tree. MST not detected in cluster tree.")
  if (is(x$hc, "hclust")) hc <- x$hc else hc <- x$hc[[which]]

  # Save default plot settings
  .pardefault <- par(no.readonly = T)

  ## Scatter plot / nodes
  if (missing(f) || is.null(f)){
    x_r2 <- list()
    if (x$d == 1){
      x_r2 <- cbind(X_n, X_n) # Mirror points
    } else if(x$d == 2){
      x_r2 <- X_n
    } else {
      ## Plot first two principle components
      pr_comp <- prcomp(as.matrix(X_n))
      x_r2 <- pr_comp$x[, 1:2]
    }
    plot(x_r2, pch = 20)
  } else {
    f(1:x$n)
  }

  ## Edges
  mst <- x[["mst"]]
  X_n <- attr(x, "X_n")
  if (is.null(h)){
    x_coords <- if (missing(f)) X_n[mst[, 1], 1] else f(mst[, 1], ...)
    y_coords <- if (missing(f)) X_n[mst[, 2], 1] else f(mst[, 2], ...)
    lines(x_coords, y_coords)
    all_pts <- if (missing(f)) X_n else f(1:nrow(X_n))
    points(all_pts, pch = 20)
  }

  # else {
  #   i <- (mst[, 3] < h)
  #
  #   cl <- cutree(hc, h = h)
  #   cl_colors <- sample(rainbow(length(unique(cl))))
  #   C_n$mst[C_n$mst[, 3] < h,]
  #   x_coords <- c(x_r2[mst[i, 1]+1, 1], x_r2[mst[i, 2]+1, 1])
  #   y_coords <- c(x_r2[mst[i, 1]+1, 2], x_r2[mst[i, 2]+1, 2])
  #   lines(x_coords, y_coords)
  #   points(x_r2, pch = 20, col = cl_colors[cl])
  # }




  par(.pardefault) # restore
}