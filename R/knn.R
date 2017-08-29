#' @title K-nearest neighbor search
#' @name knn
#' @description Performs a k-nearest neighbor search.
#' @param x a matrix-coercible data object.
#' @param k the number of neighbors to find.
#' @param r_x the reference data set, if any. See details.
#' @param bucketSize maximum size of the kd-tree leaves.
#' @param splitRule rule to split the kd-tree. One of "STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR" or "SUGGEST" (SL stands for sliding).
#' @param approx whether to allow an arbitrary small amount of error in the result for better speed.
#' @param ... unused.
#' @description "SUGGEST" uses the ANN authors preferred 'best' rule. If no reference data set (\code{rx}) is specified,
#' then the k-nearest neighbors in \code{x} are found. If a reference data set is specified, a kd tree is built using the
#' reference set, and then the k-nearest reference neighbors to each query pt are returned.
#' @references 1. Mount, David M., and Sunil Arya. "ANN: library for approximate nearest neighbour searching." (1998).
#' @export
knn <- function(x, k, r_x = NULL, bucketSize = 15L, splitRule = "suggest", approx = 0, ...){

  ## Choose the split rule for the kd tree
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  ## Check type of input for x
  if (is.null(dim(x))) stop("'knn' expects x to be a matrix-coercible object.")
  if (!is.null(r_x) && is.null(dim(r_x))) stop("'knn' expects r_x to be a matrix-coercible object.")
  x <- as.matrix(x)
  if (!storage.mode(x) %in% c("double", "integer")) stop("'clustertree' expects x to be numeric or integer only.")

  ## Future work
  # if (method == "dual tree" || missing(method)){
  #   ## Allow metric parameters to be passed via ...
  #   metric_par <- append(list(d = ncol(x)), list(...))
  #   metric_ptr <- chooseMetric(metric, metric_par)
  #   res <- dt_kNN_int(x, k = k, bucketSize = bucketSize, splitRule = splitRule, metric_ptr = metric_ptr)
  #   return(res);
  # }
  #else if (method == "single tree" && metric == "euclidean"){
  if (is.null(r_x)){
    r_x <- matrix()
    res <- kNN_int(x, k = k, type = 1, bucketSize = bucketSize, splitRule = splitRule, approx = approx, r_x = matrix(0))
  } else {
    res <- kNN_int(x, k = k, type = 1, bucketSize = bucketSize, splitRule = splitRule, approx = approx, r_x = r_x)
  }
  return(res)
}