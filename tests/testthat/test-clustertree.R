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
  expect_warning(clustertree(X_n, k = 5L, warn = T))
})

## Make sure setup code at least works
test_that("Basic functionality works", {
  expect_silent(clustertree(X_n, k = 5L))
})

## Logic check: does k = 2 w/ alpha = 1 yield the same results as single linkage
test_that("RSL with k = 2 and alpha = 1 is the same as single linkage", {
  cl <- clustertree(X_n, k = 2L, alpha = 1, estimator = "RSL")
  sl <- hclust(dist_x, method = "single")

  ## The heights should be identical
  expect_equal(cl$hc$height, sl$height, tolerance = sqrt(.Machine$double.eps))

  ## As should the cuts
  expect_equivalent(cutree(cl$hc, k = 1:(n - 1)), cutree(sl, k = 1:(n - 1)))
})

## Check euclidean distance objects produce the same results
test_that("'clustertree' gives the same results for a 'dist' objects", {
  cl <- clustertree(X_n, k = 5L, alpha = sqrt(2), estimator = "RSL")
  cl2 <- clustertree(dist(X_n, method = "euclidean"), k = 5L, alpha = sqrt(2), estimator = "RSL")

  ## The heights should be identical
  expect_equal(cl$hc$height, cl2$hc$height, tolerance = sqrt(.Machine$double.eps))

  ## As should the cuts
  expect_equivalent(cutree(cl$hc, k = 1:(n - 1)), cutree(cl2$hc, k = 1:(n - 1)))
})


## Check that KNN estimator works fine
test_that("'clustertree' works with 'KNN' estimator", {
  expect_silent(clustertree(X_n, k = 5L, alpha = sqrt(2), estimator = "KNN"))
})

## Check that mutual KNN estimator works fine
test_that("'clustertree' works with 'mutual KNN' estimator", {
  expect_silent(clustertree(X_n, k = 5L, alpha = sqrt(2), estimator = "mutual KNN"))
})

## Check that estimators which produce multiple CCs work
test_that("'clustertree' works with multiple components", {
  expect_silent(clustertree(X_n, k = 2L, alpha = sqrt(2), estimator = "mutual KNN"))
  cl <- clustertree(X_n, k = 2L, alpha = sqrt(2), estimator = "mutual KNN")
  expect_true(is.list(cl$hc)) ## expect list of hierarchies
  expect_equal(length(cl$hc), 6L)
  expect_true(all(sapply(cl$hc, function(hc) is(hc, "hclust")))) ## expect all are hierarchies
  expect_output(print(cl)) ## expect printing works
})




