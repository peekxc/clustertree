# simplified_hclust.R
# Contains a variety of S3 methods that specialize 'simplified' hclust objects, which are
# valid hierarchical 'hclust' objects where, unlike regular hclust objects, leaves correspond to >= 1
# point indices from the original data set. In this representation, branches have

plot.simplified_hclust <- function(x){

}


# cutree <- function(x) UseMethod("cutree")

#' @export
cut.simplified_hclust <- function(x, k = NULL, h = NULL){
  if (!is(x, "simplified_hclust")){ stop("'cut.simplified_hclust' expects a simplified hclust object with branch labels."); }
  if (!is(x, "hclust")){ stop("'cut.simplified_hclust' expect the simplified hclust object is a valid hclust object."); }
  # Rely on cutrees testing of inputs first
  cl <- cutree(x, k=k, h=h)
  n <- sum(sapply(ls(x$idx), function(key) length(x$idx[[key]])))
  simpl_cut <- cut_simplified_hclust(x, cl_in = cl, big_n = n)
  return(simpl_cut)
}