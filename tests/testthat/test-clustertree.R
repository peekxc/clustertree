library("clustertree")
library("testthat")

context("clustertree")

## Load stock dividend data
load(system.file("test_data/StockDividends.rdata", package = "clustertree"))
X_n <- StockDividends[, 2:12]

## Supply Variables/parameters
{ n <- nrow(X_n); dist_x <- dist(X_n, method = "euclidean") }

## Automatically detecting k ~ dlogn sets k around 39, more neighbors than records!
expect_error(clustertree::clustertree(X_n))

## RSL clustertree solution using brute-force/naive version
cl_tree_mst <- clustertree::clustertree(X_n, k = 5L)
r_k <- apply(dbscan::kNNdist(X_n, k = cl_tree_mst$k - 1), 1, max)
cl_tree_naive <- clustertree:::naive_clustertree(dist_x, r_k, cl_tree_mst$alpha, type = 0)

## Do the run time lengths of each symbol match?
expect_true(all(sapply(cl_tree_naive, function(cl) {
  mst_rle <- rle(cutree(cl_tree_mst, h = cl$r))$lengths
  bf_rle <- rle(cl$cluster)$lengths
  all(mst_rle == bf_rle)
})))

## Do the individual symbols all map to the same value, and in the correct order?
expect_true(all(sapply(cl_tree_naive, function(cl) {
  mst_vals <- rle(cutree(cl_tree_mst, h = cl$r))$values
  bf_vals <- rle(cl$cluster)$values
  all(match(mst_vals, unique(mst_vals)) == match(bf_vals, unique(bf_vals)))
})))


# which_changed <- function(cl1, cl2){
#   from_changed <- which(cl1 != cl2)
#   to_changed <- cl2[from_changed]
#   cbind(from_changed, to_changed)
# }

## Algorithm 2 Validation tests
cl_tree_knn_naive <- clustertree:::naive_clustertree(dist_x, r_k, cl_tree_mst$alpha, type = 1)

cl_tree_mst_knn <- clustertree:::primsRSL(dist_x, r_k = r_k, n = n, alpha = sqrt(2), type = 1)
cl_tree_mst_knn <- clustertree:::mstToHclust(cl_tree_mst_knn, n)

#sort(unique(sapply(cl_tree_knn_naive[-1], function(cl) cl$R))) == sort(unique(cl_tree_mst_knn$height))

# sapply(cl_tree_knn_naive, function(cl) cl$R)
expect_true(all(sapply(cl_tree_knn_naive, function(cl) {
  mst_rle <- rle(cutree(cl_tree_mst_knn, h = cl$R))$lengths
  bf_rle <- rle(cl$cluster)$lengths
  ifelse(length(mst_rle) == length(bf_rle), all(mst_rle == bf_rle), FALSE)
})))

## Do the individual symbols all map to the same value, and in the correct order?
expect_true(all(sapply(cl_tree_knn_naive, function(cl) {
  mst_vals <- rle(cutree(cl_tree_mst_knn, h = cl$R))$values
  bf_vals <- rle(cl$cluster)$values
  ifelse(length(mst_vals) == length(bf_vals),
         all(match(mst_vals, unique(mst_vals)) == match(bf_vals, unique(bf_vals))), FALSE)
})))
#
#
# data.table::rbindlist(lapply(1:15, function(i)
#   as.data.frame(which_changed(cl_tree_knn_naive[[i]]$cluster, cl_tree_knn_naive[[i+1]]$cluster))),
#   idcol = T)
#
# ## clustertree solution using MST
# cl_tree_mst_knn <- clustertree::clusterTree(dist_x = dist_x, r_k = apply(knn$dist, 1, max),
#                                             cl_tree_mst$k, alpha = cl_tree_mst$alpha,
#                                             knn_indices = knn$id[, ncol(knn$id)], type = 1)
#
# clustertree:::.mstToHclust(cl_tree_mst_knn[[1]][1:146,], n = 150)
#
# ## Do the run time lengths of each symbol match?
# expect_true(all(sapply(cl_tree_knn_naive, function(cl) {
#   mst_rle <- rle(cutree(cl_tree_mst_knn, h = cl$r))$lengths
#   bf_rle <- rle(cl$cluster)$lengths
#   ifelse(length(mst_rle) == length(bf_rle), all(mst_rle == bf_rle), FALSE)
# })))
#
#
# ## RSL clustertree solution using brute-force/naive version
# r_k <- apply(dbscan::kNNdist(X_n, k = cl_tree_mst$k - 1), 1, max)
# cl_tree_naive <- clustertree:::naive_clustertree(dist_x, r_k, cl_tree_mst$alpha)
#
#
# cl_tree_mknn_naive <- clustertree:::naive_clustertree(dist_x, r_k, cl_tree_mst$alpha, type = 2)
#
#
