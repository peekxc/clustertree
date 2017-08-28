#' @title Runt Pruning
#' @name runt_prune
#' @description Implementation of Stuetzle's 'Runt pruning' for cluster tree and generic hierarchical cluster representations.
#' @param C_n Either a 'clustertree' or a generic 'hclust' object.
#' @param runt_size Positive integer. The number of original observations a branch needs to contain to be considered "real" or "significant". See details.
#' @return If a 'clustertree' object was given, the underlying hierarchy (or hierarchies) are replaced by pruned hierarchies. If a regular
#' 'hclust' is given, a simplified (runt-pruned) version of the hierarchy is returned.
#' @description Runt pruning simplifies a given cluster tree into a new, pruned cluster tree. Algorithmically, it
#' begins by divisively traversing a given cluster tree, and at each split, looking at
#' whether a branch contains a sufficient number of observations (determined by the \code{runt_size}). If it does, the
#' branch (and its associated connected component) is considered 'significant', and is retained in the simplified tree.
#' If it does not, the entire branch and its descendents are grouped together, and are collectively recorded as either:
#' 1) A leaf, if the other branch also has less than then \code{runt_size} number of observations.
#' 2) As part of the splitting branch, if the other branch has at least \code{runt_size} number of observations.
#' @seealso cut.simplified_hclust
#' @references 1. Stuetzle, Werner. "Estimating the cluster tree of a density by analyzing the minimal spanning tree of a sample." Journal of classification 20.1 (2003): 025-047.
#' @export
runt_prune <- function(C_n, runt_size){
  if (!any(c("clustertree", "hclust") %in% class(C_n))) stop("'runt_prune' expects a 'clustertree' or 'hclust' object.")
  if (runt_size <= 0){ stop("runt size must be positive.") }
  if (runt_size == 1){ warning("Runt size must be >= 2."); return(invisible(C_n)); }

  # If it's a clustertree object, loop through the (potentially several) hierarchies and create a pruning hierarchy
  # for each one, and then return the augmented clustertree object. Otherwise, just return the pruned hclust object.
  if (is(C_n, "clustertree")){
    if (is(C_n$hc, "hclust")) { C_n$hc <- simplified_hclust(C_n$hc, runt_size) }
    else if (is(C_n$hc, "list") && length(C_n$hc) > 0) {
      C_n$hc <- lapply(C_n$hc, FUN = function(hcl) simplified_hclust(hcl, runt_size))
    }
    return(C_n)
  }
}