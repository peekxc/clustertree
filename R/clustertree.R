#' @title Compute the Empirical Cluster Tree
#' @name clustertree
#' @description More details coming soon...
#' @param x a matrix-coercible data set, or a euclidean 'dist' object.
#' @param k integer; smoothing parameter, controls the radius of the ball containing k points around each x_i.
#' @param alpha float; regularity parameter determing when to connect points. See below for details.
#' @note The default 'k' parameter will differ between if a 'dist' object is given vs. the original data set, as with
#' the dist object the original dimensionality of the data is unknown (to which the default setting of k depends on).
#' @references See KC and SD.
#' @import dbscan
#' @importFrom methods is
#' @useDynLib clustertree
#' @export
clustertree <- function(x, k = "suggest", alpha = "suggest", estimator = c("RSL", "knn", "mKnn")){
  if (is(x, "dist")){
    if (attr(x, "method") != "euclidean")
      warning("Robust Single Linkage expects euclidean distances. See ?clustertree for more details.")
    dist_x <- x
    k <- ifelse(missing(k), log(nrow(dist_x)), k)

    # original_symbol <- as.character(attr(dist_x, "call")[["x"]])
    ## Attempt to retrieve original data set and thus dimensionality
    # if (original_symbol %in% ls(parent.frame(1))){
    #   x <- as.matrix(eval(original_symbol, envir = parent.frame(1)))
    #   k <- ifelse(missing(k),  ncol(x) * log(nrow(x)), k)
    # } else { k <- ifelse(missing(k), log(nrow(x)), k) }
  } else {
    x <- as.matrix(x)
    dist_x <- dist(x, method = "euclidean")
    k <- ifelse(missing(k), ncol(x) * log(nrow(x)), k)
    if (!is.matrix(x)) stop("clustertree expects the data to be either matrix-convertible or a dist-object.")
  }
  k <- as.integer(k)
  alpha <- ifelse(missing(alpha), sqrt(2), alpha)

  ## Choose estimator
  type <- ifelse(missing(estimator), 0, pmatch(estimator, c("RSL", "knn", "mKnn")))
  if (is.na(type)){
    stop(paste0("Unknown estimator supplied. Please use one of: [", paste0(c("RSL", "knn", "mKnn"), collapse = ", "), "]"))
  }

  ## Warn about parameter settings yielding unknown results
  if (k < floor(ncol(x) * log(nrow(x))))
    warning("Existing analysis on RSL rely on alpha being at least sqrt(2) and k being at least as large as d*logn.")

  r_k <- dbscan::kNNdist(x, k = k - 1)
  res <- clusterTree(dist_x = dist_x, r_k = apply(r_k, 1, max), k = k, alpha = alpha, type = type)
  res$call <- match.call()
  res$method <- "robust single linkage"
  res$k <- k
  res$alpha <- alpha
  res
}
