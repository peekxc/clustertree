library("clustertree")
library("testthat")

context("clustertree")


dist_x <- as.matrix(dist(iris[, 1:4]))
edge_list <- t(mapply(function(i, j) c(i, j, dist_x[i, j]), row(dist_x), col(dist_x)))
mst_truth <- optrees::msTreePrim(nodes = 1:nrow(iris), arcs = edge_list)

## Test base prims algorithm
iris_mst <- clustertree:::primsMST(dist(iris[, 1:4]))


sum(mst_truth$tree.arcs[, 3])
testthat::expect_equal(sum(mst_truth$tree.arcs[, 3]), sum(iris_mst[, 3]))


clustertree:::dtb()