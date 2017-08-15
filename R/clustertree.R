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
clustertree <- function(x, d = ncol(x), k = "suggest", alpha = "suggest",
                        estimator = c("RSL", "knn", "mutual knn"),
                        warn = FALSE){
  if (is(x, "dist")){
    d <- ifelse(missing(d), NULL, d)
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
  warn_message <- "Existing clustertree analysis relies on alpha being at least sqrt(2) and k being at least as large as d*log(n)."
  if (k < ceiling(ncol(x) * log(nrow(x))) && warn) warning(warn_message)

  r_k <- dbscan::kNNdist(x, k = k - 1)
  hc <- clusterTree(dist_x = dist_x, r_k = apply(r_k, 1, max), k = k, alpha = alpha, type = type)
  hc$call <- match.call()
  hc$method <- possible_estimators[type+1]
  res <- structure(list(hc = hc, k = k, d = d, alpha = alpha), class = "clustertree")
  return(res)
}

#' @export
print.clustertree <- function(C_n){
  type <- pmatch(C_n$hc$method, c("RSL", "knn", "mutual knn"))
  est_type <- c("Robust Single Linkage", "KNN graph", "Mutual KNN graph")[type+1]
  writeLines(c(
    paste0("clustertree object estimated using: ", est_type),
    sprintf("Parameters: k = %d, alpha = %.3f, dim = %d", C_n$k, C_n$alpha, C_n$d)
  ))
}

