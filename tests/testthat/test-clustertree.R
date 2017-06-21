library("clustertree")
library("testthat")

context("clustertree")

## Load stock dividend data
load(system.file("test_data/StockDividends.rdata", package = "clustertree"))
X_n <- StockDividends[, 2:12]

## Supply Variables/parameters
{ n <- nrow(X_n); dist_x <- dist(X_n, method = "euclidean") }

## Automatically detecting k ~ dlogn sets k around 39, more neighbors than records!
test_that("Parameter detection works", {
    expect_error(clustertree::clustertree(X_n))
    expect_warning(clustertree::clustertree(X_n, k = 5L, warn_parameter_settings = T))
})

## Make sure setup code at least works
test_that("Basic functionality works", {
  testthat::expect_silent(clustertree::clustertree(X_n, k = 5L))
  testthat::expect_silent(apply(dbscan::kNNdist(X_n, k = 5L - 1), 1, max))
})

## RSL clustertree solution using brute-force/naive version
test_that("RSL naive solution matches MST", {
  mst_rsl <- clustertree::clustertree(X_n, k = 5L)
  r_k <- apply(dbscan::kNNdist(X_n, k = mst_rsl$k - 1), 1, max)
  rsl_naive <- clustertree:::naive_clustertree(dist_x, r_k, mst_rsl$alpha, type = 0)

  ## Do the run time lengths of each symbol match?
  expect_true(all(sapply(rsl_naive, function(cl) {
    mst_rle <- rle(cutree(mst_rsl, h = cl$R))$lengths
    bf_rle <- rle(cl$cluster)$lengths
    all(mst_rle == bf_rle)
  })))

  ## Do the individual symbols all map to the same value, and in the correct order?
  expect_true(all(sapply(rsl_naive, function(cl) {
    mst_vals <- rle(cutree(mst_rsl, h = cl$R))$values
    bf_vals <- rle(cl$cluster)$values
    all(match(mst_vals, unique(mst_vals)) == match(bf_vals, unique(bf_vals)))
  })))
})

## Algorithm 2 Validation tests
test_that("kNN naive solution matches MST", {
  mst_knn <- clustertree::clustertree(X_n, k = 5L, estimator = "knn")
  r_k <- apply(dbscan::kNNdist(X_n, k = mst_knn$k - 1), 1, max)
  knn_naive <- clustertree:::naive_clustertree(dist_x, r_k, mst_knn$alpha, type = 1)

  expect_true(all(sapply(knn_naive, function(cl) {
    mst_rle <- rle(cutree(mst_knn, h = cl$R))$lengths
    bf_rle <- rle(cl$cluster)$lengths
    ifelse(length(mst_rle) == length(bf_rle), all(mst_rle == bf_rle), FALSE)
  })))

  ## Do the individual symbols all map to the same value, and in the correct order?
  expect_true(all(sapply(knn_naive, function(cl) {
    mst_vals <- rle(cutree(mst_knn, h = cl$R))$values
    bf_vals <- rle(cl$cluster)$values
    ifelse(length(mst_vals) == length(bf_vals),
           all(match(mst_vals, unique(mst_vals)) == match(bf_vals, unique(bf_vals))), FALSE)
  })))
})

test_that("mutual kNN naive solution matches MST", {
  mst_mknn <- clustertree::clustertree(X_n, k = 5L, estimator = "knn")
  r_k <- apply(dbscan::kNNdist(X_n, k = mst_mknn$k - 1), 1, max)
  mknn_naive <- clustertree:::naive_clustertree(dist_x, r_k, mst_mknn$alpha, type = 1)

  expect_true(all(sapply(mknn_naive, function(cl) {
    mst_rle <- rle(cutree(mst_mknn, h = cl$R))$lengths
    bf_rle <- rle(cl$cluster)$lengths
    ifelse(length(mst_rle) == length(bf_rle), all(mst_rle == bf_rle), FALSE)
  })))

  ## Do the individual symbols all map to the same value, and in the correct order?
  expect_true(all(sapply(mknn_naive, function(cl) {
    mst_vals <- rle(cutree(mst_mknn, h = cl$R))$values
    bf_vals <- rle(cl$cluster)$values
    ifelse(length(mst_vals) == length(bf_vals),
           all(match(mst_vals, unique(mst_vals)) == match(bf_vals, unique(bf_vals))), FALSE)
  })))
})
