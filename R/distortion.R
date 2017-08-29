#' Merge distortion
#' @param hc1 hclust object; the first hierarchy
#' @param hc2 hclust object; the second hierarchy
#' @param correspondence the correspondence function to use.
#' @description Implements the merge distortion metric described in [1]. Currently only supports the identity correspondence function.
#' @references 1. Eldridge, Justin, Mikhail Belkin, and Yusu Wang. "Beyond hartigan consistency: Merge distortion metric for hierarchical clustering." Conference on Learning Theory. 2015.
#' @export
distortion <- function(hc1, hc2, correspondence = "identity"){
  if (!is(hc1, "hclust")) stop("'distortion' metric expects two valid hclust objects!")
  if (!is(hc2, "hclust")) stop("'distortion' metric expects two valid hclust objects!")
  if (nrow(hc1$merge) != nrow(hc2$merge)){
    stop("Distortion metric expects the hierarchies to be built from the same data sets!")
  }
  m1 <- mergeHeight(hc1)
  m2 <- mergeHeight(hc2)
  return(max(abs(m1 - m2)))
}