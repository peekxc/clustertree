#' @title Plot a given cluster tree
#' @name plot.clustertree
#' @description More details coming soon...
#' @param x a 'clustertree' object.
#' @references See KC and SD.
#' @importFrom methods is
#' @useDynLib clustertree
#' @export
plot.clustertree <- function(x){

}




# TODO: Check for dimensionality, store the mst in the clustertree object
spanplot <- function(X_n, C_n, h = NULL){
  if (!any(c("clustertree") %in% class(C_n))) stop("spanplot expects a 'clustertree' object")
  hc <- ifelse(is(C_n, "clustertree"), C_n$hc, C_n)
  if (is.null(h)){
    plot(X_n, pch = 20)
    for (i in 1:nrow(C_n$height)){

    }
    cutree(C_n, h = 2)
  } else {
    plot(X_n)
  }

}