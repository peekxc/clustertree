#' @title K-nearest neighbor search
#' @name clustertree
#' @description Performs a k-nearest neighbor search.
#' @param x a matrix-coercible data object.
#' @param k the number of neighbors to find.
#' @param r_x the reference data set, if any. See details.
#' @param metric the metric to use.
#' @param method whether to use a dual tree traversal (any supported metric) or a single tree traversal (euclidean only).
#' @param bucketSize maximum size of the kd-tree leaves.
#' @param splitRule rule to split the kd-tree. One of "STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR" or "SUGGEST" (SL stands for sliding).
#' "SUGGEST" uses the ANN authors preferred 'best' rule.
#' @references 1. Mount, David M., and Sunil Arya. "ANN: library for approximate nearest neighbour searching." (1998).
#' @export
knn <- function(x, k, r_x = NULL, metric = "euclidean", method = c( "single tree"),
                bucketSize = 15L, splitRule = "suggest",
                ...){
  if (method == "single tree" && metric != "euclidean") stop("Single tree method only supports euclidean metric.")

  ## Choose the split rule for the kd tree
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  ## Check type of input for x
  if (is.null(dim(x))) stop("'clustertree' expects x to be a matrix-coercible object.")
  x <- as.matrix(x)
  if (!storage.mode(x) %in% c("double", "integer")) stop("'clustertree' expects x to be numeric or integer only.")

  if (method == "dual tree" || missing(method)){
    ## Allow metric parameters to be passed via ...
    metric_par <- append(list(d = ncol(x)), list(...))
    metric_ptr <- chooseMetric(metric, metric_par)
    res <- dt_kNN_int(x, k = k, bucketSize = bucketSize, splitRule = splitRule, metric_ptr = metric_ptr)
    return(res);
  }
  else if (method == "single tree" && metric == "euclidean"){
    res <- kNN_int(x, k = k, type = 1, bucketSize = bucketSize, splitRule = splitRule, approx = 0)
  }
  else { stop("Unknown method specified.") }
}