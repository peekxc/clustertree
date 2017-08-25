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
  type <- ifelse(missing(estimator), 0, pmatch(toupper(estimator), toupper(possible_estimators)) - 1L)
  if (is.na(type)){
    stop(paste0("Unknown estimator supplied. Please use one of: [", paste(possible_estimators, collapse = ", "), "]"))
  }

  ## Warn about parameter settings yielding unknown results
  warn_message <- "Existing clustertree analysis relies on alpha being at least sqrt(2) and k being at least as large as d*log(n)."
  if (k < ceiling(d * log(nrow(x))) && warn) warning(warn_message)

  ## Call the cluster tree function
  hc <- clusterTree_int(x = x, k = k, alpha = alpha, type = type)

  ## If it's a list, dealing with a minimum spanning forest
  if (is(hc, "list")){
    for(i in 1:length(hc)){
      if (is(hc[[i]], "hclust")){
        hc[[i]]$call <- match.call()
        hc[[i]]$method <- possible_estimators[type+1]
      }
      # Else do nothing - singleton
    }
  } else {
  ## Otherwise it's a fully connected tree
    if (is(hc, "hclust")){
      hc$call <- match.call()
      hc$method <- possible_estimators[type+1]
    }
  }
  res <- structure(list(hc = hc, k = k, d = d, alpha = alpha), class = "clustertree", X_n = x)
  return(hc)
}

#' print.clustertree
#' @method print clustertree
#' @export
print.clustertree <- function(C_n){
  type <- pmatch(toupper(C_n$hc$method), toupper(c("RSL", "knn", "mutual knn")))
  est_type <- c("Robust Single Linkage", "KNN graph", "Mutual KNN graph")[type+1]
  writeLines(c(
    paste0("Cluster tree object of: ", nrow(attr(C_n, "X_n")), " objects."),
    paste0("Estimator used: ", est_type),
    sprintf("Parameters: k = %d, alpha = %.4f, dim = %d", C_n$k, C_n$alpha, C_n$d)
  ))
}

#' @title Plot a given cluster tree
#' @name plot.clustertree
#' @description More details coming soon...
#' @param x a 'clustertree' object.
#' @references See KC and SD.
#' @importFrom methods is
#' @useDynLib clustertree
#' @export
plot.clustertree <- function(C_n, type = c("both", "dendrogram", "span tree")){
  X_n <- attr(C_n, "X_n")
  # dev.interactive()
  if (is(C_n$hc, "hclust")){
    layout(matrix(c(1, 2), nrow = 1, ncol = 2))
    plot(C_n$hc, hang = -1)
    spanplot(X_n, C_n)
  } else if (is(C_n$hc, "list") && length(C_n$hc) > 1){
    layout(matrix(c(1, 1:length(C_n$hc) +1), nrow = 1, ncol = length(C_n$hc)+1),
           widths = c(rep(0.5/length(C_n$hc), length(C_n$hc)), 0.5))
    for (hcl in C_n$hc){ plot(hcl, hang = -1) }
    spanplot(X_n, C_n)
  }
}


## The metrics supported by the various dual tree extensions
.supported_metrics <- c("euclidean", "manhattan", "maximum", "minkowski")

## The various splitting routines for the ANN kd trees
.ANNsplitRule <- c("STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR", "SUGGEST")

.onUnload <- function (libpath) {
  library.dynam.unload("clustertree", libpath)
}
