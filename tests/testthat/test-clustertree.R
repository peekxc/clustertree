library("clustertree")
library("testthat")

data("iris")
context("clustertree")

## Sample from iris data set
# X_n <- iris[1:100, 1:4]
# n <- nrow(X_n)
# k <- 5L


library("dbscan")
data(moons)
X_n <- moons
n <- nrow(X_n)
k <- 5L

r_k <- apply(dbscan::kNNdist(X_n, k = k - 1), 1, max)


test_indices <- clustertree:::checkKruskals(dist(X_n), r_k = r_k, k = 5L)

cl_tree_r <- clustertree_ex(X_n, k = 5L, alpha = sqrt(2))

## which steps were new branches created
sapply(clustertree[-1], function(branch) branch$iter)
test_indices$i
x_2 <- iris[1:50, 1:4]
what1 <- clustertree_ex(x_2, k = 5L, alpha = sqrt(2))
what2 <- clustertree_ex(x_2, k = 5L, alpha = 1)
## Brute-force R solution can be computed using igraph and FNN packages
# Created with: cl_tree_r <- clustertree_ex(X_n, k = 5L, alpha = sqrt(2))
load(system.file("test_data/iris_cl_tree_r.rdata", package = "clustertree"))

## Rcpp solution using MST
cl_tree_cpp <- clustertree(X_n, k = 5L, alpha = sqrt(2))

## Verify brute-force type solution matches C++ RSL solution
expect_true(all(sapply(cl_tree_r, function(cltree_cut) all(cutree(cl_tree_cpp, h = cltree_cut$radius) == cltree_cut$cluster))))

test_indices <- do.call(rbind, lapply(order(as.vector(dist(X_n)))-1, function(i) what(i, 100)))


cl_tree_r[[2]]$radius
r_k <- apply(dbscan::kNNdist(X_n, k = k - 1), 1, max)
test_indices <- clustertree:::checkKruskals(dist(X_n), r_k = r_k, k = 5L)

test_indices$i
sapply(cl_tree_r, function(branch) branch$iter)

test <- matrix(0, ncol= 10, nrow=10)
indices <- do.call(rbind, lapply(0:(choose(10, 2) - 1), function(i) what(i, n = 10)))[, 2:1]
c_i <- 1
for (i in 1:nrow(indices)){
  test[indices[i, 1]+1, indices[i, 2]+1] <- c_i
  c_i <- c_i + 1
}

what <- function(k, n){
  i =  n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
  j =  k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
  c(i, j)
}

distX <- as.vector(dist(X_n))
test_indices <- clustertree:::checkKruskals(dist(X_n), r_k = r_k, k = 5L)

truth <- do.call(rbind, lapply(order(as.vector(dist(X_n))) - 1, function(i) what(i, 100)))
#do.call(rbind, lapply(as.vector(test[lower.tri(test)]) - 1, function(i) what(i, 10)))


