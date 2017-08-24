library("clustertree")
library("testthat")

context("clustertree")

## Load stock dividend data
load(system.file("test_data/StockDividends.rdata", package = "clustertree"))
X_n <- StockDividends[, 2:12]


sl <- hclust(dist(as.matrix(X_n)), method = "single")

k <- 2
sl_simple <- clustertree:::simplified_hclust(sl, 3)

k_cuts <- 1:nrow(sl$merge)
possible_cl <- cutree(sl, k = k_cuts)

filter_noise <- function(cl){
  cl_freq <- table(cl)
  cl_ids <- as.character(cl)
  unname(ifelse(cl_freq[cl_ids] >= 2, as.integer(names(cl_freq[cl_ids])), 0))
}

apply(apply(possible_cl, 2, filter_noise), 2, function(cl) cl[sl$order])