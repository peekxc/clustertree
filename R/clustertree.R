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
clustertree <- function(x, k = "suggest", alpha = "suggest",
                        estimator = c("RSL", "knn", "mutual knn"),
                        warn = FALSE){
  if (is.null(dim(x))) stop("'clustertree' expects x to be a matrix-coercible object.")
  x <- as.matrix(x)
  if (!storage.mode(x) %in% c("double", "integer")) stop("'clustertree' expects x to be numeric or integer only.")
  d <- ncol(x)
  k <- as.integer(ifelse(missing(k), d * log(nrow(x)), k)) # dim should work now
  alpha <- ifelse(missing(alpha), sqrt(2), alpha)

  ## Choose estimator
  possible_estimators <- c("RSL", "knn", "mutual knn")
  type <- ifelse(missing(estimator), 0, pmatch(estimator, possible_estimators) - 1L)
  if (is.na(type)){
    stop(paste0("Unknown estimator supplied. Please use one of: [", paste(possible_estimators, collapse = ", "), "]"))
  }

  ## Warn about parameter settings yielding unknown results
  warn_message <- "Existing clustertree analysis relies on alpha being at least sqrt(2) and k being at least as large as d*log(n)."
  if (k < ceiling(d * log(nrow(x))) && warn) warning(warn_message)

  ## Call the cluster tree function
  # r_k <- knn(x, k = k - 1)
  hc <- clusterTree_int(x = x, k = k, alpha = alpha, type = type)
  hc$call <- match.call()
  hc$method <- possible_estimators[type+1]
  res <- structure(list(hc = hc, k = k, d = d, alpha = alpha), class = "clustertree")
  return(res)
}

print.clustertree <- function(C_n){
  type <- pmatch(C_n$hc$method, c("RSL", "knn", "mutual knn"))
  est_type <- c("Robust Single Linkage", "KNN graph", "Mutual KNN graph")[type+1]
  writeLines(c(
    paste0("clustertree object estimated using: ", est_type),
    sprintf("Parameters: k = %d, alpha = %.3f, dim = %d", C_n$k, C_n$alpha, C_n$d)
  ))
}

#' @export
plot.clustertree <- function(C_n, type = c("both", "dendrogram", "span tree")){

}




## The metrics supported by the various dual tree extensions
.supported_metrics <- c("euclidean", "manhattan", "maximum", "minkowski")

## The various splitting routines for the ANN kd trees
.ANNsplitRule <- c("STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR", "SUGGEST")

.onUnload <- function (libpath) {
  library.dynam.unload("clustertree", libpath)
}
