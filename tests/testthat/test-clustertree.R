library("clustertree")
library("testthat")

data("iris")
context("clustertree")

## Sample from iris data set
X_n <- iris[, 1:4]

## Supply Variables/parameters
{ n <- nrow(X_n); dist_x <- dist(X_n, method = "euclidean") }

## clustertree solution using MST
cl_tree_mst <- clustertree::clustertree(X_n)

## RSL clustertree solution using brute-force/naive version
knn <- dbscan::kNN(X_n, k = cl_tree_mst$k - 1)
cl_tree_naive <- clustertree:::naive_clustertree(dist_x, apply(knn$dist, 1, max), cl_tree_mst$alpha, type = 0)

## Do the run time lengths of each symbol match?
expect_true(all(sapply(cl_tree_naive, function(cl) {
  mst_rle <- rle(cutree(cl_tree_mst, h = cl$r))$lengths
  bf_rle <- rle(cl$cluster)$lengths
  all(mst_rle == bf_rle)
})))

## Do the individual symbols all map to the same value, an in the correct order?
expect_true(all(sapply(cl_tree_naive, function(cl) {
  mst_vals <- rle(cutree(cl_tree_mst, h = cl$r))$values
  bf_vals <- rle(cl$cluster)$values
  all(match(mst_vals, unique(mst_vals)) == match(bf_vals, unique(bf_vals)))
})))

# ## Algorithm 2 Validation tests
# cl_tree_knn_naive <- clustertree:::naive_clustertree(dist_x, apply(knn$dist, 1, max), cl_tree_mst$alpha, type = 1)
#
# ## clustertree solution using MST
# cl_tree_mst_knn <- clustertree::clusterTree(dist_x = dist_x, r_k = apply(knn$dist, 1, max),
#                                             cl_tree_mst$k, alpha = cl_tree_mst$alpha,
#                                             knn_indices = knn$id[, ncol(knn$id)], type = 1)
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


