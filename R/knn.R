#' @title K-nearest neighbor search
#' @name clustertree
#' @description More details coming soon...
#' @param x a matrix-coercible data object.
#' @param k the number of neighbors to find.
#' @param r_x the reference data set, if any. See details.
#' @param metric the metric to use.
#' @param method whether to use a dual tree traversal (any supported metric) or a single tree traversal (euclidean only).
#' @param bucketSize maximum size of the kd-tree leaves.
#' @param splitRule rule to split the kd-tree. One of "STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR" or "SUGGEST" (SL stands for sliding).
#' "SUGGEST" uses the ANN authors preferred 'best' rule.
#' @export
knn <- function(x, k, r_x = NULL, metric = "euclidean", method = c("dual tree", "single tree"),
                bucketSize = 15L,
                splitRule = "suggest"){
  # pmatch(metric, c(""))
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")
  if (method == "dual tree" || missing("dual tree")){
    # Run the kNN
    res <- dt_kNN_int(X_n, k = k, bucketSize = bucketSize, splitRule = 5L, metric_ptr = test)
    return(res)
  } else if (method == "single tree" && metric == "euclidean"){
    # Run the kNN
    metric_ptr
    res <- dt_kNN_int(X_n, k = k, bucketSize = bucketSize, splitRule = 5L, metric_ptr = test)
    return(res)
  }

}