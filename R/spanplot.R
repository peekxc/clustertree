#' @title Plot a spanning tree
#' @name spanplot
#' @description Allows for visualizing a clustertree object as a spanning tree.
#' @param X_n a matrix to plot as a scatter plot.
#' @param C_n a 'clustertree' object
#' @param h an optional height to cut the tree at. Colors the corresponding connected components.
#' @importFrom methods is
#' @export
spanplot <- function(X_n, C_n, h = NULL){
  if (!any(c("clustertree") %in% class(C_n))) stop("spanplot expects a 'clustertree' object")
  hc <- C_n$hc
  if (is.null(hc[["mst"]])) stop("Cannot plot spanning tree. MST not detected in cluster tree.")

  ## Scatter plot / nodes
  x_r2 <- list()
  if (C_n$d == 1){
    x_r2 <- cbind(X_n, X_n) # Mirror points
  } else if(C_n$d == 2){
    x_r2 <- X_n
  } else {
    ## Plot first two principle components
    pr_comp <- stats::prcomp(as.matrix(X_n))
    x_r2 <- pr_comp$x[, 1:2]
  }
  plot(x_r2, pch = 20)

  ## Edges
  if (is.null(h)){
    for (i in 1:nrow(hc$mst)){
      x_coords <- c(x_r2[hc$mst[i, 1]+1, 1], x_r2[hc$mst[i, 2]+1, 1])
      y_coords <- c(x_r2[hc$mst[i, 1]+1, 2], x_r2[hc$mst[i, 2]+1, 2])
      lines(x_coords, y_coords)
    }
    points(x_r2, pch = 20)
  } else {
    cl <- cutree(hc, h = h)
    cl_colors <- sample(rainbow(length(unique(cl))))
    for (i in 1:nrow(hc$mst)){
      if (hc$height[i] >= h) { break }
      x_coords <- c(x_r2[hc$mst[i, 1]+1, 1], x_r2[hc$mst[i, 2]+1, 1])
      y_coords <- c(x_r2[hc$mst[i, 1]+1, 2], x_r2[hc$mst[i, 2]+1, 2])
      lines(x_coords, y_coords)
    }
    points(x_r2, pch = 20, col = cl_colors[cl])
  }
}