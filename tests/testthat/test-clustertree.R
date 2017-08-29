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
  expect_true(all(cl$hc$height == sl$height))

  ## As should the cuts
  expect_true(all(cutree(cl$hc, k = 1:(n - 1)) ==  cutree(sl, k = 1:(n - 1))))
})
