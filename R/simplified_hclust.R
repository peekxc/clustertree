# simplified_hclust.R
# Contains a variety of S3 methods that specialize 'simplified' hclust objects, which are
# valid hierarchical 'hclust' objects where, unlike regular hclust objects, leaves correspond to >= 1
# point indices from the original data set. In this representation, branches have

#' @importFrom stats hclust
#' @import graphics
#' @export
plot.simplified_hclust <- function(x, ...){
  if (nrow(x$merge) == 1){
    ## Hclust doesn't handle very very simple hierarchies very well, so convert to dendrogram and plot that instead
    plot(as.dendrogram(x))
  }else {
    # Swap order of classes to allow hclust plotting
    class(x) <- c("hclust", "simplified_hclust")
    plot(x)
  }
}

#' @export
cut.simplified_hclust <- function(x, ...){
  if (!is(x, "simplified_hclust")){ stop("'cut.simplified_hclust' expects a simplified hclust object with branch labels."); }
  if (!is(x, "hclust")){ stop("'cut.simplified_hclust' expect the simplified hclust object is a valid hclust object."); }
  config <- list(...)
  h <- NULL
  k <- NULL
  if (!is.na(pmatch("h", names(config)))){ h <- config[["h"]] }
  if (!is.na(pmatch("k", names(config)))){ k <- config[["k"]] }
  if (is.null(k) && is.null(h)){ stop("Please specify a height or cut value to cut the hierarchy at."); }
  if (!is.null(k) && !is.null(h)){ stop("You can only specify the number of clusters to extract (k) or a height (h), not both."); }

  ## Rely on cutrees testing of inputs for further checking
  cl <- cutree(x, k=k, h=h)

  ## Make the cut using the result from cutree
  n <- sum(sapply(ls(x$idx), function(key) length(x$idx[[key]])))
  simpl_cut <- cut_simplified_hclust(x, cl_in = cl, big_n = n)
  return(simpl_cut)
}