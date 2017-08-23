#' @title Compute the Empirical Cluster Tree
#' @name clustertree
#' @description Implements various estimators for estimating the empirical cluster tree.
#' @param x a matrix-coercible data set, or a euclidean 'dist' object.
#' @param k integer; smoothing parameter, controls the radius of the ball containing k points around each x_i.
#' @param alpha float; regularity parameter that prevents excessive chaining. Can be safely kept on the default. See below for details.
#' @note The default 'k' parameter will differ between if a 'dist' object is given vs. the original data set, as with
#' the dist object the original dimensionality of the data is unknown (to which the default setting of k depends on).
#' Chaudhuri et. al recommend keeping alpha as low as possible, subject to being no lower than sqrt(2).
#' @references See Chaudhuri, Kamalika, and Sanjoy Dasgupta. "Rates of convergence for the cluster tree." Advances in Neural Information Processing Systems. 2010.
#' @import dbscan
#' @importFrom methods is
#' @useDynLib clustertree
#' @export
clustertree <- function(x, k = "suggest", alpha = "suggest", estimator = c("RSL", "knn", "mutual knn"),
                        warn_parameter_settings = FALSE){
  if (is(x, "dist")){
    if (attr(x, "method") != "euclidean" && warn_parameter_settings)
      warning("Robust Single Linkage expects euclidean distances. See ?clustertree for more details.")
    dist_x <- x
    k <- ifelse(missing(k), log(nrow(dist_x)), k)
  } else {
    x <- as.matrix(x)
    dist_x <- dist(x, method = "euclidean")
    k <- ifelse(missing(k), ncol(x) * log(nrow(x)), k)
    if (!is.matrix(x)) stop("clustertree expects the data to be either matrix-convertible or a dist-object.")
  }
  k <- as.integer(k)
  alpha <- ifelse(missing(alpha), sqrt(2), alpha)

  ## Choose estimator
  possible_estimators <- c("RSL", "knn", "mutual knn")
  type <- ifelse(missing(estimator), 0, pmatch(estimator, possible_estimators) - 1L)
  if (is.na(type)){
    stop(paste0("Unknown estimator supplied. Please use one of: [", paste(possible_estimators, collapse = ", "), "]"))
  }

  ## Warn about parameter settings yielding unknown results
  if (k < floor(ncol(x) * log(nrow(x))) && warn_parameter_settings)
    warning("Existing analysis on RSL rely on alpha being at least sqrt(2) and k being at least as large as d*logn.")

  r_k <- dbscan::kNNdist(x, k = k - 1)
  res <- clusterTree(dist_x = dist_x, r_k = apply(r_k, 1, max), k = k, alpha = alpha, type = type)
  res$call <- match.call()
  res$method <- "robust single linkage"
  res$k <- k
  res$alpha <- alpha
  res
}

## The metrics supported by the various dual tree extensions
.supported_metrics <- c("euclidean", "manhattan", "maximum", "minkowski")

## The various splitting routines for the ANN kd trees
.ANNsplitRule <- c("STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR", "SUGGEST")

.onUnload <- function (libpath) {
  library.dynam.unload("clustertree", libpath)
}
