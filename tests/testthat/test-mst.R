library("clustertree")
library("testthat")

context("mst")

X_n <- as.matrix(iris[, 1:4])
sl <- hclust(dist(X_n), method = "single")

## The euclidean mininum spanning tree is equivalent to single linkage. Thus, the
## sum of the edge weights from both should be the same
euc_mst <- mst(dist(X_n))
expect_true(sum(sl$height) == sum(euc_mst[, 3]))


