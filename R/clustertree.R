#' @title Compute the Empirical Cluster Tree
#' @name clustertree
#' @description Implements various algorithms for estimating the cluster tree.
#' @param x a matrix-coercible data set, or a euclidean 'dist' object.
#' @param k integer; smoothing parameter controlling how 'dense' an object should be before being considered for
#' linkage into the cluster tree hierarchy. See below for details.
#' @param alpha float; regularity parameter used to prevent 'excessive chaining'. Can be safely kept on the default. See below for details.
#' @param estimator Which estimator to choose from. Default is robust single linkage, see below for details.
#' @param warn boolean; Whether to warn the user regarding the usage of the parameter settings. Defaults to false. See below for details.
#' @note The default 'k' parameter will differ between if a 'dist' object is given vs. the original data set, as with
#' the dist object the original dimensionality of the data is unknown (to which the default setting of k depends on).
#' Chaudhuri et. al [1] recommend keeping alpha as low as possible, but recommended no lower than sqrt(2). If not estimator
#' is specified, Robust Single Linkage (RSL) is used.
#'
#' The other two estimators are from Algorithm 2 in [2]. Note that in the case where the KNN or Mutual KNN graph estimators are used,
#' the resulting 'tree' may actually be a disconnected forest. To handle this each connected component is made into its
#' own hierarchy ('hclust' object), discarding singletons (if any). Because \code{hclust} objects require a representation
#' where leaf indices start with 1, the indices of the original data set are stored in the \code{idx} member of each
#' hierarchy.
#' @references
#' 1. Chaudhuri, Kamalika, and Sanjoy Dasgupta. "Rates of convergence for the cluster tree." Advances in Neural Information Processing Systems. 2010.
#' 2. Chaudhuri, Kamalika, et al. "Consistent procedures for cluster tree estimation and pruning." IEEE Transactions on Information Theory 60.12 (2014): 7900-7912.
#' @seealso hclust
#' @importFrom methods is
#' @importFrom Rcpp evalCpp
#' @useDynLib clustertree
#' @export
clustertree <- function(x, k = "suggest", alpha = "suggest",
                        estimator = c("RSL", "knn", "mutual knn"),
                        warn = FALSE){
  if (is.null(dim(x)) && !is(x, "dist")) stop("'clustertree' expects x to be a matrix-coercible object.")
  if (!is(x, "dist")) x <- as.matrix(x) ## Coerce to matrix
  if (!storage.mode(x) %in% c("double", "integer")) stop("'clustertree' expects x to be numeric or integer only.")
  x_dist <- NULL
  n <- 0
  if (is(x, "dist")){
    d <- 1L
    n <- attr(x, "Size")
    if (warn && x$method != "euclidean"){ warning("Current estimators have only been studied under the euclidean metric.") }
    x_dist <- x
  } else {
    d <- ncol(x)
    n <- nrow(x)
    if (requireNamespace("parallelDist", quietly = TRUE)){
      x_dist <- parallelDist::parallelDist(x, method = "euclidean") ## prefer parallel version of dist if possible
    } else { x_dist <- dist(x, method = "euclidean") }
  }
  k <- as.integer(ifelse(missing(k), d * log(n), k)) # dim should work now
  alpha <- ifelse(missing(alpha), sqrt(2), alpha)

  r_k <- NULL
  if (is(x, "dist")){
    r_k <- knn_dist2(x_dist, k - 1)
  } else {
    r_k <- knn(x, k = k - 1, bucketSize = k * log(n), splitRule = "SUGGEST")
    r_k <- apply(r_k$dist, 1, max)
  }

  ## Make sure k is less than n
  if (k > n){
    if (warn) warning("Data set has either too few records or to high of a dimensionality to use recommended default parameters.")
    k <- max(2, ceiling(log(n)))
  }

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
  ## Note: knn returns the k nearest, exclusively. Clustertree relies on k to be inclusive, so increase by 1.
  st <- primsRSL(x_dist, r_k = r_k, n = n, alpha = alpha, type = type)

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
    # hclust_info <- mstToHclust(st[, 1:2], st[, 3])
    hclust_info <- hclustMergeOrder(st, order(st[, 3]))
    hclust_info$call <- match.call()
    hclust_info$method <- possible_estimators[type+1]
  }
  mst <- st[order(st[, 3]),]
  mst <- st[!st[, 3] %in% c(.Machine$double.xmax, NA),]
  res <- structure(list(hc = hclust_info, k = k, d = d, n = n, alpha = alpha, mst = mst), class = "clustertree", X_n = x)
  return(res)
}

#' print.clustertree
#' @description Prints various parameters used to compute the cluster tree.
#' @param x a 'clustertree' object.
#' @param ... unused.
#' @method print clustertree
#' @export
print.clustertree <- function(x, ...){
  method <- ifelse(is(x$hc, "list"), x$hc[[1]]$method, x$hc$method)
  type <- pmatch(toupper(method), c("RSL", "KNN", "MUTUAL KNN"))
  est_type <- c("Robust Single Linkage", "KNN graph", "Mutual KNN graph")[type]
  writeLines(c(
    paste0("Cluster tree object of: ", x$n, " objects."),
    paste0("Estimator used: ", est_type),
    sprintf("Parameters: k = %d, alpha = %.4f, dim = %d", x$k, x$alpha, x$d)
  ))
}

#' @title Plot a given cluster tree
#' @name plot.clustertree
#' @description More details coming soon...
#' @param x a 'clustertree' object.
#' @param which which hierarchy to plot. If RSL was used or k was sufficiently small, this may be only 1 hierarchy.
#' @param h height to cut the tree(s) at in the visualization. If null, all of the connected components are visualized.
#' @param type whether to plot the clustertree where the height is computed using the radius at which clusters form (root at the top),
#' or whether to use the corresponding empirical density (root at the bottom).
#' @param ... arguments to pass to the plotting function. Otherwise unused.
#' @description This plotting function plots the hierarchy specified by the \code{which} argument that makes up part of
#' the cluster tree object on the left half of plotting device and the spanning tree of the points within that hierarchy on
#' the right.
#' @seealso \code{spanplot}
#' @import graphics
#' @export
plot.clustertree <- function(x, which = 1, h = NULL, type=c("radius", "density"), ...){
  X_n <- attr(x, "X_n")
  # dev.interactive()
  .pardefault <- par(no.readonly = T)

  hc <- NULL
  if (is(x$hc, "hclust")){ hc <- x$hc }
  else { hc <- x$hc[[which]] }
  if (!is(hc, "hclust")) stop("No valid cluster tree hierarchy found.")

  est_type <- pmatch(toupper(hc$method), c("RSL", "KNN", "MUTUAL KNN"))
  estimator <- "Cluster tree"
  if (!is.na(est_type)){
    estimator <- c("Robust Single Linkage", "KNN Graph", "Mutual KNN Graph")[est_type]
  }

  ## Plot the dendrogram with the appropriate title
  if (missing(type) || type == "radius"){
    plot(hc, ylab = expression(italic(r)), main = estimator, ...)
  } else if (type == "density"){
    r <- hc$height
    v_d <- vol_nSphere(x$d)
    lambda <- x$k / (x$n * v_d * r^x$d)
    x$hc$height <- lambda ## Replace height with lambda
    plot(as.dendrogram(x$hc), ylab = expression(lambda), main = estimator,
         ylim = rev(range(lambda)), ...)
  } else {
    stop("Unknown type supplied.")
  }


  # if (is(x$hc, "hclust")){
  #   layout(matrix(c(1, 2), nrow = 1, ncol = 2))
  #   plot(x$hc, hang = -1)
  #   if (is.null(h)){
  #     spanplot(X_n, x)
  #   } else {
  #     spanplot(X_n, x, h = h)
  #   }
  # } else if (is(x$hc, "list") && length(x$hc) > 0){
  #   n_hier <- length(x$hc)
  #   if (which > n_hier){
  #     stop(sprintf("Cannot plot hierarchy (%d), there are only %d total hierarchies!", which, n_hier))
  #   }
  #   layout(matrix(c(1, 2), nrow = 1, ncol = 2))
  #   plot(x$hc[[which]], hang = -1)
  #   if (missing(f) || is.null(f)){
  #     spanplot(X_n, x)
  #   } else {
  #     f(X_n, x, ...)
  #   }
  # }
  par(.pardefault)
}

# Generic way to apply an operation to every valid 'hclust' obejct in the clustertree.
# applyCT <- function(C_n, expr, f){
#   if (is(C_n$hc, "hclust")){
#
#   } else if (is(C_n$hc, "list") && length(C_n$hc) > 0){
#
#   }
# }

# TODO
# identify.clustertree <- function(C_n, which = 1, h = NULL, f="scatter"){
#
# }

## The metrics supported by the various dual tree extensions
.supported_metrics <- c("euclidean", "manhattan", "maximum", "minkowski")

## The various splitting routines for the ANN kd trees
.ANNsplitRule <- c("STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR", "SUGGEST")

.onUnload <- function (libpath) {
  library.dynam.unload("clustertree", libpath)
}
