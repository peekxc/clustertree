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
  x <- as.matrix(x) ## Coerce to matrix
  if (!storage.mode(x) %in% c("double", "integer")) stop("'clustertree' expects x to be numeric or integer only.")
  x_dist <- NULL
  if (is(x, "dist")){
    d <- 1L
    if (warn && x$method != "euclidean"){ warning("Current estimators have only been studied under the euclidean metric.") }
    x_dist <- x
  } else {
    d <- ncol(x)
    x_dist <- dist(x, method = "euclidean")
  }
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
  if (warn && k < ceiling(d * log(nrow(x)))) warning(warn_message)

  ## Call the cluster tree augmented MST (using prioritized prims with delayed connections)
  r_k <- knn(x, k = k, bucketSize = k * log(nrow(x)), splitRule = "SUGGEST")
  r_k <- apply(r_k$dist, 1, max)
  st <- primsRSL(x_dist, r_k = r_k, n = nrow(x), alpha = alpha, type = type)

  ## If it's from estimators 1 or 2, it could be a minimum spanning forest
  hclust_info <- list()
  if (type > 0){
    CCs <- mstToCC(st[, 1:2], st[, 3])
    hclust_info <- vector(mode = "list", length = sum(table(CCs) >= 2))
    i <- 1
    for(cc in unique(CCs)){
      ids <- which(CCs == cc) - 1 # Point ids in this connected component (0-based)
      if (length(ids) >= 2){
        mst_idx <- union(which(st[, 1] %in% ids), which(st[, 2] %in% ids)) # MST indices that make this CC
        sub_st <- st[mst_idx, ] # spanning component subset
        sub_st <- sub_st[order(sub_st[, 3]),] # ordered by height
        sub_st <- sub_st[!sub_st[, 3] %in% c(.Machine$double.xmax, NA),] # discard special forest linkage flags
        hc <- mstToHclust(sub_st[, 1:2], sub_st[, 3])
        hc$call <- match.call()
        hc$method <- possible_estimators[type+1]
        hc$idx <- ids + 1 # original point indices from x used to create this hierarchy
        hclust_info[[i]] <- structure(hc, class="hclust") # Enforce class
        i <- i + 1
      }
      # Else do nothing - need at least two points to make a cluster
    }
  } else {
  ## RSL creates a fully connected hierarchy, convert to an hclust object
    hclust_info <- mstToHclust(st[, 1:2], st[, 3])
    hclust_info$call <- match.call()
    hclust_info$method <- possible_estimators[type+1]
  }
  mst <- st[order(st[, 3]),]
  mst <- st[!st[, 3] %in% c(.Machine$double.xmax, NA),]
  res <- structure(list(hc = hclust_info, k = k, d = d, alpha = alpha, mst = mst), class = "clustertree", X_n = x)
  return(res)
}

#' @name print.clustertree
#' @method print clustertree
#' @export
print.clustertree <- function(C_n){
  method <- ifelse(is(C_n$hc, "list"), what$hc[[1]]$method, C_n$hc$method)
  type <- pmatch(toupper(method), c("RSL", "KNN", "MUTUAL KNN"))
  est_type <- c("Robust Single Linkage", "KNN graph", "Mutual KNN graph")[type]
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
  } else if (is(C_n$hc, "list") && length(C_n$hc) > 0){
    n_hier <- length(C_n$hc)
    layout(matrix(c(1, 1:n_hier +1), nrow = 1, ncol = n_hier+1),
           widths = c(rep(0.5/n_hier,n_hier), 0.5))
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
