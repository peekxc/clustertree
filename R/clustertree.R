#' @title Compute the Empirical Cluster Tree
#' @name clustertree
#' @description More details coming soon...
#' @param x a matrix-coercible data set, or a euclidean 'dist' object.
#' @param k integer; smoothing parameter, controls the radius of the ball containing k points around each x_i.
#' @param alpha float; regularity parameter determing when to connect points. See below for details.
#' @note The default 'k' parameter will differ between if a 'dist' object is given vs. the original data set, as with
#' the dist object the original dimensionality of the data is unknown (to which the default setting of k depends on).
#' @references See KC and SD.
#' @export
clustertree <- function(x, k = ncol(x) * log(nrow(x)), alpha = sqrt(2)){
  if (is(x, "dist")){
    if (attr(x, "method") != "euclidean")
      warning("Robust Single Linkage expects euclidean distances. See ?clustertree for more details.")
    dist_x <- x
    x <- as.matrix(eval(attributes(dist_x)$call[["x"]], envir = parent.env(environment())))
    # k <- ifelse(missing(k), log(attr(dist_x, "Size")), k)
  } else {
    x <- as.matrix(x)
    dist_x <- dist(x, method = "euclidean")
    if (!is.matrix(x)) stop("clustertree expects the data to be either matrix-convertible or a dist-object.")
  }
  k <- as.integer(k)
  r_k <- dbscan::kNNdist(x, k = k - 1)
  res <- clusterTree(x = dist_x, r_k = r_k[, k - 1], k = k, alpha = alpha)
  res$call <- match.call()
  res$method <- "robust single linkage"
  res
}
