library("clustertree")
library("testthat")

context("knn")

## Load stock dividend data
load(system.file("test_data/StockDividends.rdata", package = "clustertree"))
X_n <- StockDividends[, 2:12]

## Make sure knn doesn't error
expect_silent(knn(X_n, k = 10L))

## Calculate truth
all_dist <- as.matrix(dist(as.matrix(X_n), method = "euclidean"))
knn_truth <- list(
  id = t(apply(all_dist, 1, order))[, 2:11],
  dist = t(apply(all_dist, 1, sort))[, 2:11],
  k = 10
)

## Run the k-nearest neighbor + validate with dist
knn_res <- knn(X_n, k = 10L)
expect_equivalent(knn_res$id, knn_truth$id)
expect_true(all.equal(knn_res$dist, knn_truth$dist, tolerance = sqrt(.Machine$double.eps), check.attributes = FALSE))


