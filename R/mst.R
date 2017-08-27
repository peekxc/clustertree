#' Minimum spanning tree
#' @param X_n a 'dist' object or matrix of edge weights encoding an adjacency matrix.
#' @details Given a 'dist' object encoding the effective 'edge weights' between nodes, calculate the
#' minimum spanning tree. Nonexistent edges can be represented as \code{NA} or \code{.Machine$double.xmax}. Only
#' supports symmetric edge weights.
#' @return a matrix containing the 'from' node indices (column 1), the 'to' node indices (column 2), and the weight of
#' the edge between then (column 3).
#' @export
mst <- function(X_n){
  if (!storage.mode(X_n) %in% c("integer", "double")){ stop("mst expects an integer or numeric type matrix of weights.") }
  if (is.null(dim(X_n)) && !is(X_n, "dist")) stop("mst expects a 'dist' object or a weighted adjacency matrix.")
  if (is(X_n, "matrix")){
    X_n <- as.dist(X_n) # coerce a 'dist' object
  }
  res <- primsMST(X_n)
  return(res)
}