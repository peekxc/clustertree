#' Minimum spanning tree
#' @description Given a data set, calculate the euclidean minimum spanning tree. Equivalent to single linkage, but
#' avoids calculating the distance matrix.
#' @export
mst <- function(X_n, metric = "euclidean"){
  x <- as.matrix(X_n)
  metric_ptr <- chooseMetric(metric, config = list(d = ncol(X_n)))
  res <- dtb_int(x, metric_ptr)
  return(res)
}