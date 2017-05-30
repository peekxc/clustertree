library("clustertree")
library("testthat")

data("iris")
context("clustertree")

## Sample from iris data set
X_n <- iris[1:100, 1:4]
n <- nrow(X_n)

## Brute-force R solution can be computed using igraph and FNN packages
# cl_tree_r <- clustertree_ex(X_n, k = 5L, alpha = sqrt(2))
load(system.file("test_data/iris_cl_tree_r.rdata", package = "clustertree"))

## Rcpp solution using MST
cl_tree_cpp <- clustertree(X_n, k = 5L, alpha = sqrt(2))

## Verify brute-force type solution matches C++ RSL solution
expect_true(all(sapply(cl_tree_r, function(cltree_cut) all(cutree(cl_tree_cpp, h = cltree_cut$radius) == cltree_cut$cluster))))