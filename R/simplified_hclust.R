# simplified_hclust.R
# Contains a variety of S3 methods that specialize 'simplified' hclust objects, which are
# valid hierarchical 'hclust' objects where, unlike regular hclust objects, leaves correspond to >= 1
# point indices from the original data set. In this representation,

plot.simplified_hclust <- function(x){

}


# cutree <- function(x) UseMethod("cutree")

#' @export
cut.simplified_hclust <- function(x, k = NULL, h = NULL){

  # Rely on cutrees testing of inputs first
  cl <- cutree(x, k=k, h=h)

  n <- sum(sapply(ls(what$idx), function(key) length(what$idx[[key]])))
  cl_out <- vector(mode = "integer", length = n)
  idx
  print("hello")
}
