// R_dt.cpp
// Function to group up internal entry-point / API calls for all of the dual tree algorithms

#include "DualTreeSearch.h"
#include "DualTreeKNN.h"
#include "DualTreeBoruvka.h"

// [[Rcpp::export]]
void testTrees(const NumericMatrix& x, const int k, const int bucketSize = 15){
  // Copy data over to ANN point array
  ANNpointArray queryPts = matrixToANNpointArray(x);

  // Test creating a regular tree
  ANNkd_tree* kdTree1 = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5);

  // Test creating a dual-tree
  NODE_INFO node_info; // container to store various properties related to tree construction
  ANNkd_tree* kdTree2 = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5, node_info);
  // Rcout << "Seconds to create dual tree: " << (t1 - t0) / 1000000.0L << std::endl;

  delete kdTree1;
  delete kdTree2;
}

// [[Rcpp::export]]
List testClusterTree(const NumericMatrix& x, const int k, const int bucketSize = 15, const double alpha = 1.0){ //1.41421356237
  // Copy data over to ANN point array
  ANNpointArray queryPts = matrixToANNpointArray(x);
  const int n = x.nrow();

  // Create the kd tree
  NODE_INFO node_info; // container to store various properties related to tree construction
  ANNkd_tree* kdTreeR = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5, node_info);

  // Get the KNN neighbors
  L_2* l2_norm = new L_2(queryPts, queryPts, x.ncol());
  DualTreeKNN<L_2> dts = DualTreeKNN<L_2>(kdTreeR, kdTreeR, node_info, node_info, *l2_norm, k); // -1 to be inclusive of point itself
  dts.DFS(); // Depth-First Traversal solves the KNN

  // Get the KNN distance of the ball containing k neighbors
  std::vector<ANNdist> rk_dist = std::vector<ANNdist>();
  rk_dist.reserve(n);
  for (int i = 0; i < n; ++i){
    rk_dist.push_back(dts.knn_pq[i]->max_key());
  }

  // Create the robust single linkage metric space
  RSL* rsl_metric = new RSL(queryPts, queryPts, x.ncol(), alpha, rk_dist);

  // Rprintf("RSL Dist: 0 --> 1 : %f, alpha: %f, r_k(0): %f, r_k(1): %f\n", rsl_metric->distance(0, 1), rsl_metric->alpha, rsl_metric->r_k[0], rsl_metric->r_k[1]);
  std::chrono::time_point<std::chrono::high_resolution_clock> start_total = std::chrono::high_resolution_clock::now();
  DualTreeBoruvka<RSL> dtb = DualTreeBoruvka<RSL>(kdTreeR, kdTreeR, node_info, node_info, *rsl_metric);
  NumericMatrix mst = dtb.MST();
  std::chrono::time_point<std::chrono::high_resolution_clock> finish_total = std::chrono::high_resolution_clock::now();

  Rprintf("Number of pts compared: %d\n", dtb.n_comparisons);
  Rprintf("Number of saved pts accessed: %d\n", dtb.n_accesses);
  Rprintf("Number of pruned recursions: %d\n", dtb.n_pruned);

  // Print the benchmarks
  double total_time = std::chrono::duration<double, std::milli>(finish_total - start_total).count();
  Rprintf("Total time: %f ms \n", total_time);
  double recorded_percent = 0.0;
  for (int i = 0; i < benchmarks.size(); ++i){
    // std::chrono::duration_cast<std::chrono::seconds>(benchmarks.at(i));
    Rprintf("Benchmark %d, %-100s: %f ms (%.2f%%)\n", i, benchmark_descrip.at(i).c_str(), benchmarks.at(i), (benchmarks.at(i)/total_time)*100 );
    recorded_percent += (benchmarks.at(i)/total_time)*100;
  }
  Rprintf("Percent accounted for (%.2f%%)\n", recorded_percent);


  // Make the return result
  List res = List::create(_["mst"] = mst, _["tree_properties"] = dtb.getTreeProperties(), _["knn"] = wrap(rk_dist));
  delete kdTreeR;
  return(res);
}

// [[Rcpp::export]]
List testKNN(const NumericMatrix& x, const int k, const int bucketSize = 15, const int metric_choice = 0){
  // Copy data over to ANN point array
  ANNpointArray queryPts = matrixToANNpointArray(x);

  // Create the kd tree
  NODE_INFO node_info = NODE_INFO(); // container to store various properties related to tree construction
  ANNkd_tree* kdTreeR = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5, node_info);

  // Rcout << "node info size: " << node_info.size() << std::endl;
  // for (NODE_INFO::iterator ni = node_info.begin(); ni != node_info.end(); ++ni){
  //   Rcout << ni->id << std::endl;
  // }
  // return(List::create());

  // Create the Dual Tree KNN object, pointing to the new node info
  List res;
  if (metric_choice == L2){
    L_2* l2_norm = new L_2(queryPts, queryPts, x.ncol());
    DualTreeKNN<L_2> dts = DualTreeKNN<L_2>(kdTreeR, kdTreeR, node_info, node_info, *l2_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == L1){
    L_1* l1_norm = new L_1(queryPts, queryPts, x.ncol());
    DualTreeKNN<L_1> dts = DualTreeKNN<L_1>(kdTreeR, kdTreeR, node_info, node_info, *l1_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == LINF){
    L_inf* linf_norm = new L_inf(queryPts, queryPts, x.ncol());
    DualTreeKNN<L_inf> dts = DualTreeKNN<L_inf>(kdTreeR, kdTreeR, node_info, node_info, *linf_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == LP){
    L_p* lp_norm = new L_p(queryPts, queryPts, x.ncol(), 2);
    DualTreeKNN<L_p> dts = DualTreeKNN<L_p>(kdTreeR, kdTreeR, node_info, node_info, *lp_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == RSL_M){
    // Get the KNN neighbors
    L_2* l2_norm = new L_2(queryPts, queryPts, x.ncol());
    DualTreeKNN<L_2> dts = DualTreeKNN<L_2>(kdTreeR, kdTreeR, node_info, node_info, *l2_norm, k - 1); // -1 to be inclusive of point itself
    dts.DFS(); // Depth-First Traversal solves the KNN

    // Get the KNN distance of the ball containing k neighbors
    std::vector<ANNdist> rk_dist = std::vector<ANNdist>();
    rk_dist.reserve(x.nrow());
    for (int i = 0; i < x.nrow(); ++i){
      rk_dist.push_back(dts.knn_pq[i]->max_key());
    }

    // Create the robust single linkage metric space
    RSL* rsl_metric = new RSL(queryPts, queryPts, x.ncol(), 1, rk_dist);
    DualTreeKNN<RSL> dts2 = DualTreeKNN<RSL>(kdTreeR, kdTreeR, node_info, node_info, *rsl_metric, k);
    dts2.DFS(); // Depth-First Traversal solves the KNN
    res = dts2.getKnnResults();
    res["tree_properties"] = dts2.getTreeProperties();
    delete kdTreeR;
  }

  return(res);
}

// [[Rcpp::export]]
List testKNN_ref(const NumericMatrix& qx, const NumericMatrix& rx, const int k, const int bucketSize = 15, const int metric_choice = 0){
  if (qx.ncol() != rx.ncol()){ stop("Query point dimension != reference point dimension."); }

  // Copy data over to ANN point array
  ANNpointArray queryPts = matrixToANNpointArray(qx);
  ANNpointArray refPts = matrixToANNpointArray(rx);

  // Create the query kd tree
  NODE_INFO node_info_q = NODE_INFO(); // container to store various properties related to tree construction
  ANNkd_tree* kdTreeQ = new ANNkd_tree(queryPts, qx.nrow(), rx.ncol(), bucketSize, (ANNsplitRule) 5, node_info_q);

  NODE_INFO node_info_r = NODE_INFO(); // container to store various properties related to tree construction
  ANNkd_tree* kdTreeR = new ANNkd_tree(refPts, rx.nrow(), rx.ncol(), bucketSize, (ANNsplitRule) 5, node_info_r);

  // for (NODE_INFO::iterator ni = node_info_q.begin(); ni != node_info_q.end(); ++ni){
  //   Rprintf("Q Node id: %d\n", ni->id);
  // }
  // Rcout << "Q Number of nodes: " << node_info_q.size() << std::endl;
  // for (NODE_INFO::iterator ni = node_info_r.begin(); ni != node_info_r.end(); ++ni){
  //   Rprintf("R Node id: %d\n", ni->id);
  // }
  // Rcout << "R Number of nodes: " << node_info_r.size() << std::endl;

  // Create the Dual Tree KNN object, pointing to the new node info
  List res;
  if (metric_choice == L2){
    L_2* l2_norm = new L_2(queryPts, refPts, qx.ncol());
    DualTreeKNN<L_2> dts = DualTreeKNN<L_2>(kdTreeQ, kdTreeR, node_info_q, node_info_r, *l2_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == L1){
    L_1* l1_norm = new L_1(queryPts, refPts, qx.ncol());
    DualTreeKNN<L_1> dts = DualTreeKNN<L_1>(kdTreeQ, kdTreeR, node_info_q, node_info_r, *l1_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == LINF){
    L_inf* linf_norm = new L_inf(queryPts, refPts, qx.ncol());
    DualTreeKNN<L_inf> dts = DualTreeKNN<L_inf>(kdTreeQ, kdTreeR, node_info_q, node_info_r, *linf_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == LP){
    L_p* lp_norm = new L_p(queryPts, refPts, qx.ncol(), 2);
    DualTreeKNN<L_p> dts = DualTreeKNN<L_p>(kdTreeQ, kdTreeR, node_info_q, node_info_r, *lp_norm, k);
    dts.DFS(); // Depth-First Traversal solves the KNN
    res = dts.getKnnResults();
    res["tree_properties"] = dts.getTreeProperties();
    delete kdTreeR;
  } else if (metric_choice == RSL_M){
    // Get the KNN neighbors
    L_2* l2_norm = new L_2(queryPts, refPts, qx.ncol());
    DualTreeKNN<L_2> dts = DualTreeKNN<L_2>(kdTreeQ, kdTreeR, node_info_q, node_info_r, *l2_norm, k - 1); // -1 to be inclusive of point itself
    dts.DFS(); // Depth-First Traversal solves the KNN

    // Get the KNN distance of the ball containing k neighbors
    std::vector<ANNdist> rk_dist = std::vector<ANNdist>();
    rk_dist.reserve(qx.nrow());
    for (int i = 0; i < qx.nrow(); ++i){
      rk_dist.push_back(dts.knn_pq[i]->max_key());
    }

    // Create the robust single linkage metric space
    RSL* rsl_metric = new RSL(queryPts, refPts, qx.ncol(), 1, rk_dist);
    DualTreeKNN<RSL> dts2 = DualTreeKNN<RSL>(kdTreeQ, kdTreeR, node_info_q, node_info_r, *rsl_metric, k);
    dts2.DFS(); // Depth-First Traversal solves the KNN
    res = dts2.getKnnResults();
    res["tree_properties"] = dts2.getTreeProperties();
    delete kdTreeR;
  }

  return(res);
}


// [[Rcpp::export]]
List testMST(const NumericMatrix& x, const int bucketSize = 15){
  // Copy data over to ANN point array
  ANNpointArray queryPts = matrixToANNpointArray(x);

  // Create the dual tree
  NODE_INFO node_info; // container to store various properties related to tree construction
  ANNkd_tree* kdTreeR = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5, node_info);

  // Create the Dual Tree Boruvka object, pointing to the new node info
  L_2 l2_norm = L_2(queryPts, queryPts, x.ncol());
  DualTreeBoruvka<L_2> dts = DualTreeBoruvka<L_2>(kdTreeR, kdTreeR, node_info, node_info, l2_norm);

  // Calculate the MST
  NumericMatrix mst = dts.MST();

  Rprintf("Number of pts compared: %d\n", dts.n_comparisons);
  Rprintf("Number of saved pts accessed: %d\n", dts.n_accesses);
  Rprintf("Number of pruned recursions: %d\n", dts.n_pruned);

  // Print the benchmarks
  for (int i = 0; i < benchmarks.size(); ++i){
    // std::chrono::duration_cast<std::chrono::seconds>(benchmarks.at(i));
    Rprintf("Benchmark %d: %f ms\n", i, benchmarks.at(i));
  }

  // Get result and return
  List res = List::create(_["mst"] = mst,
                          _["tree_properties"] = dts.getTreeProperties());
  // List res = List::create(_["mst"] = mst);
  return(res);
}

/*** R
# x <- as.matrix(iris[, 1:4])
set.seed(1234)
n <- 150
x <- rbind(cbind(rnorm(n/2), rnorm(n/2)), cbind(rnorm(n/2, mean = 5), rnorm(n/2, mean = 5)))
info <- clustertree:::testKNN(x, k = 5, bucketSize = 2L)
info_db <- dbscan::kNN(x, k = 5)
all(info_db$id == info$idx)
all(abs(info_db$dist^2 - info$dist) < sqrt(.Machine$double.eps))

for (i in 1:1000){
  set.seed(i)
  n <- 15
  n2 <- 30
  x <- rbind(cbind(rnorm(n/2), rnorm(n/2)), cbind(rnorm(n/2, mean = 5), rnorm(n/2, mean = 5)))
  x2 <- rbind(cbind(rnorm(n2/2), rnorm(n2/2)), cbind(rnorm(n2/2, mean = 5), rnorm(n2/2, mean = 5)))
  random_k <- sample(1:7, size = 1)
  bck_size <- sample(1:n, size = 1)
  res <- clustertree:::testKNN_ref(x, x2, k = random_k, bucketSize = bck_size)
  qr_res <- FNN::get.knnx(data = x2, query = x, k = random_k)
  if (!all(qr_res$nn.index == res$idx)){ break; }
  if (!all(abs(qr_res$nn.dist^2 - res$dist) < sqrt(.Machine$double.eps))) { break; }
}

plotTreeProp(res$tree_properties)



res <- clustertree:::testKNN_ref(x, x2, k = 5, bucketSize = 2L)


set.seed(1234)
n <- 1000
x <- rbind(cbind(rnorm(n/2), rnorm(n/2)), cbind(rnorm(n/2, mean = 5), rnorm(n/2, mean = 5)))
microbenchmark::microbenchmark({ info <- clustertree:::testKNN(x, k = 800, bucketSize = 15L) }, times = 15L)
microbenchmark::microbenchmark({ info_db <- dbscan::kNN(x, k = 800, bucketSize = 15L, sort = FALSE) }, times = 15L)
all(info_db$id == info$idx)
all(abs(info_db$dist^2 - info$dist) < sqrt(.Machine$double.eps))


n <- 1000
x <- rbind(cbind(rnorm(n/2), rnorm(n/2)), cbind(rnorm(n/2, mean = 5), rnorm(n/2, mean = 5)))
res <- clustertree:::ClusterTree_int(x, k = 5L)
clustertree:::mstToHclust(res$mst[, 1:2], res$mst[, 3])

max(c(dist(x[1:2,])^2, core_dist[1]^2, core_dist[2]^2))

test_cltree <- clustertree:::testClusterTree(x, k = 1, alpha = 1)
all(abs(test_cltree$knn - info_db$dist[, 1]^2) < sqrt(.Machine$double.eps))
sum(test_cltree$mst[, 3]) == sum(hclust(dist(x), method = "single")$height^2)


true_mst <- cbind(hclust(dist(x)^2, method = "single")$merge, hclust(dist(x)^2, method = "single")$height)

info <- clustertree:::testDTS(x, k = nrow(x)-1, bucketSize = 2L, metric_choice = 4)
as.dist(cbind(0, info$dist))

target_k <- nrow(x) - 1
alpha <- 1
core_dist <- dbscan::kNNdist(x, k = 1)[, 1]
mrd_D <- dbscan:::mrd(dist(x, method = "euclidean")^2/alpha, cd = core_dist)
attr(mrd_D, "Size") <- nrow(x)
class(mrd_D) <- "dist"
as.matrix(mrd_D)


dbscan:::prims( nrow(x))
core_dist <- dbscan::kNNdist(x, k = 1)[, 1]
wut <- clustertree:::testClusterTree(x, k = 1, bucketSize = 5L, alpha = 1)
core_dist == wut$mst[, 3]
# core_dist == wut$knn
# sort(core_dist) == sort(sqrt(wut$knn))
sl_test <- hclust(dist(x), method = "single")
sum(sqrt(wut$mst[, 3]))
sum(sl_test$height)



sum(wut$mst[, 3])


k_ <- 13
knn_truth <- dbscan::kNN(dist(x, method = "maximum"), k = k_)

## Just computing raw KNN results from distance matrix
dist_matrix <- as.matrix(dist(x, method = "maximum"))
dist_truth <- list(id = t(apply(apply(dist_matrix, 1, order), 2, function(idx) idx[2:(k_+1)])),
              dist = t(apply(dist_matrix, 1, function(dist_x) sort(dist_x)[2:(k_+1)])))

all.equal(knn_truth$id, dist_truth$id, check.attributes = FALSE)
all.equal(knn_truth$dist, dist_truth$dist, check.attributes = FALSE)

info <- clustertree:::testDTS(x, k = 13, bucketSize = 1L, metric_choice = 2)
all.equal(info$idx, dist_truth$id, check.attributes = FALSE)
all.equal(info$dist, dist_truth$dist, check.attributes = FALSE)

plot(x, xlim = range(x[, 1]) + c(-1, 1), ylim = range(x[, 2]) + c(-1, 1), asp = 1)
plotTreeProp(info$tree_properties)
text(x[, 1], x[, 2], pos =3, labels = 1:nrow(x))

test_cltree <- clustertree:::testClusterTree(x, k = 1, bucketSize = 5L, alpha = 1)
dist(x)^2

test_mst <- clustertree:::testMST(x, bucketSize = 10)
wut <- dbscan:::prims(dist(x)^2, nrow(x))
wut_x <- hclust(dist(x), method = "single")
## These should all be the same
sum(test_mst$mst[, 3])
sum(wut[, 3])
sum(wut_x$height^2)

plot(x, xlim = range(x[, 1]) + c(-1, 1), ylim = range(x[, 2]) + c(-1, 1), asp = 1)
plotTreeProp(test_mst$tree_properties)
text(x[, 1], x[, 2], pos =3, labels = 1:nrow(x))

## Stress testing
n <- 15
for (i in 1:1000){
  set.seed(i)
  n_ <- as.integer(n/2)
  x <- rbind(cbind(rnorm(n_), rnorm(n_)),
             cbind(rnorm(n_, mean = 5), rnorm(n_, mean = 5)))
  bck_sz <- as.integer(runif(1, min = 1, max = as.integer(n_/2)))
  test_sq_dist <- sum(clustertree:::testMST(x, bucketSize = 5L)$mst[, 3])
  true_sq_dist <- sum(dbscan:::prims(dist(x)^2, nrow(x))[, 3])
  if (abs(true_sq_dist - test_sq_dist) > sqrt(.Machine$double.eps)){ break; }
}
sum(test_mst$mst[, 3])
sum(wut_x$height^2)



test_mst <- clustertree:::testMST(x, bucketSize = 5L)
plot(x, xlim = range(x[, 1]) + c(-1, 1), ylim = range(x[, 2]) + c(-1, 1), asp = 1)
plotTreeProp(test_mst$tree_properties)
text(x[, 1], x[, 2], pos =3, labels = 1:nrow(x))
wut_x <- hclust(dist(x), method = "single")
# true_mst <- dbscan:::prims(dist(x)^2, nrow(x))

## Stress testing
n <- 100

naive_cltree <- function(x, k){
  core_dist <- dbscan::kNNdist(x, k = k)[, k]
  base_dist <- dbscan:::mrd(dm = (dist(x)^2)/sqrt(2), cd = core_dist^2)
  class(base_dist) <- "dist"
  attr(base_dist, "Size") <- nrow(x)
  info <- clustertree:::primsMST(base_dist)
  return(info)
}

dt_cltree <- function(x, k){
  test_cltree <- clustertree:::testClusterTree(x, k = k, alpha = sqrt(2))
  return(test_cltree)
}

dt_cltree2 <- function(x, k){
  test_cltree <- clustertree:::Cl(x, k = k, alpha = sqrt(2))
  return(test_cltree)
}

n <- 10000
n_ <- as.integer(n/2)
x <- rbind(cbind(rnorm(n_), rnorm(n_)),
           cbind(rnorm(n_, mean = 5), rnorm(n_, mean = 5)))
random_k <- sample(1:(n-1), size = 1)
peakRAM(naive_cltree(x, random_k))
peakRAM(dt_cltree(x, random_k))


true_mst <- clustertree:::primsMST(dist(x)^2)

mst_lines <- lapply(1:nrow(true_mst), function(i){ rbind(x[true_mst[i, 1], 1:2], x[true_mst[i, 2], 1:2]) })
for (mst_line in mst_lines){
  lines(x = mst_line[, 1], y = mst_line[, 2])
}



## Stress testing
n <- 100
for (i in 1:1000){
  set.seed(i)
  n_ <- as.integer(n/2)
  x <- rbind(cbind(rnorm(n_), rnorm(n_)),
             cbind(rnorm(n_, mean = 5), rnorm(n_, mean = 5)))
  bck_sz <- as.integer(runif(1, min = 1, max = as.integer(n_/2)))
  info <- clustertree:::testKNN(x, k = as.integer(n/2), bucketSize = bck_sz) ## 11,19 should be pruned
  info_db <- dbscan::kNN(x, k = as.integer(n/2))
  res1 <- all(info_db$id == info$idx)
  res2 <- all(abs(info_db$dist^2 - info$dist) < sqrt(.Machine$double.eps))
  if (res1 != TRUE || res2 != TRUE){ break; }
}




rbind(unname(info_db$id)[1,], info$idx[1,])
rbind(unname(info_db$dist^2)[1,], info$dist[1,])

info <- clustertree:::testDTS(x, k = floor(nrow(x)/5), bucketSize = bck_sz)
info_db <- dbscan::kNN(x, k = floor(nrow(x)/5))
all(info_db$id == info$idx)
all((info_db$dist^2 - info$dist) < sqrt(.Machine$double.eps))

printf <- function(...) invisible(print(sprintf(...)))
wrong_idx <- which(apply(info_db$id == info$idx, 1, function(tf) any(!tf)))
printf("Incorrect: %d := [%s]\n", wrong_idx, paste0(info$idx[wrong_idx,], collapse = ", "))
printf("Correct  : %d := [%s]\n", wrong_idx, paste0(unname(info_db$id[wrong_idx,]), collapse = ", "))

#pts: \(q:\d+, r:9

n <- 40
for (i in 1:1000){
  set.seed(i)
  n_ <- as.integer(n/2)
  x <- rbind(cbind(rnorm(n_), rnorm(n_)),
             cbind(rnorm(n_, mean = 5), rnorm(n_, mean = 5)))
  random_k <- sample(1:n_, size = 1)
  core_dist <- dbscan::kNNdist(x, k = random_k)[, random_k]
  base_dist <- dbscan:::mrd(dm = (dist(x)^2)/sqrt(2), cd = core_dist^2)
  class(base_dist) <- "dist"
  attr(base_dist, "Size") <- nrow(x)
  info <- clustertree:::primsMST(base_dist)
  test_cltree <- clustertree:::testClusterTree(x, k = random_k, alpha = sqrt(2))
  if (abs(sum(info[, 3]) - sum(test_cltree$mst[, 3])) > sqrt(.Machine$double.eps)) { break; }
}


clustertree:::testClusterTree(x, k = random_k, alpha = sqrt(2))


plot(x, xlim = range(x[, 1]) + c(-1, 1), ylim = range(x[, 2]) + c(-1, 1), asp = 1)
plotTreeProp(info$tree_properties)
text(x[, 1], x[, 2], pos =3, labels = 1:nrow(x))

idx <- 1
res <- sapply(info$tree_properties, function(node) {
  is_in <- ((x[idx,1] >= node$lb[1]) && (x[idx,1] <= node$ub[1])) && ((x[idx,2] >= node$lb[2]) && (x[idx, 2] <= node$ub[2]))
  if (is_in){
    node$id
  } else { NULL }
})
# 1: (2|19|16|14|1|17|0)
# 17: (28|29|30|31|32|26|0|24)
# DFS: (q:16 [l:17, r:20], r:29 [l:30, r:35])


points(x = x[wrong_idx,1], y = x[wrong_idx,2], pch = 20, cex = 1.25, col = "green")

which(sapply(info$tree_properties, function(node) (17) %in% node$id))
info$tree_properties[[13]]$id

which(sapply(info$tree_properties, function(node) any((c(41, 34)-1) %in% node$idx)))
info$tree_properties[[15]]$id

info <- clustertree:::testDTS(x, k = 2, bucketSize = 3) ## 11,19 should be pruned
info_db <- dbscan::kNN(x, k = 2)
all(info_db$id == info$idx)
all((info_db$dist^2 - info$dist) < sqrt(.Machine$double.eps))

## Testing MST on small data set
info <- clustertree:::testMST(x, bucketSize = 1L)
(sum(dbscan:::prims(dist(x), nrow(x))[, 3]) - sum(sqrt(info$mst[, 3]))) < sqrt(.Machine$double.eps)

info <- clustertree:::testMST(x, bucketSize = 5L)
(sum(dbscan:::prims(dist(x), nrow(x))[, 3]) - sum(sqrt(info$mst[, 3]))) < sqrt(.Machine$double.eps)

## Testing on large data set
n <- 1000
x2 <- rbind(cbind(rnorm(n), rnorm(n), rnorm(n)),
            cbind(rnorm(n, mean = 5), rnorm(n, mean = 5), rnorm(n, mean = 5)))
all(dim(unique(x2)) == dim(x2))

## Testing validity
info <- clustertree:::testDTS(x2, k = 25, bucketSize = 25L)
## 25, 25 := 1623337 comparisons, 440 pruned
## 25, 25 := 1408497 comparisons, 2709 pruned, 458216 accesses
## 25, 25 (noref_basecase) := 1566747 comparisons, 2312 pruned, 582107 accesses
info_db <- dbscan::kNN(x2, k = 25, bucketSize = 50)
all(info_db$id == info$idx)
all((info_db$dist^2 - info$dist) < sqrt(.Machine$double.eps))

## Benchmarking
# 2165410
# 2215427
microbenchmark::microbenchmark({ info <- clustertree:::testDTS(x2, k = 25, bucketSize = 30L) }, times = 15L)



## Testing MST on larger data set
info <- clustertree:::testMST(x2, bucketSize = 50L)
wut <- dbscan:::prims(dist(x2)^2, nrow(x2))
sum(wut[, 3])
sum(info$mst[, 3])
abs(sum(wut[, 3]) - sum(sqrt(info$mst[, 3]))) < sqrt(.Machine$double.eps)

## Uniform
n <- 1000L
x3 <- cbind(runif(n), runif(n))
info <- clustertree:::testDTS(x3, k = 50, bucketSize = 15L) ## 11,19 should be pruned
info_db <- dbscan::kNN(x3, k = 50, bucketSize = 15L)
all(info_db$id == info$idx)
all((info_db$dist^2 - info$dist) < sqrt(.Machine$double.eps))

plotTreeProp <- function(tree_prop, col = "red"){
  for (node in tree_prop){
    lines(x = c(node$lb[1], node$ub[1]), y = c(node$lb[2], node$lb[2]))
    lines(x = c(node$ub[1], node$ub[1]), y = c(node$lb[2], node$ub[2]))
    lines(x = c(node$ub[1], node$lb[1]), y = c(node$ub[2], node$ub[2]))
    lines(x = c(node$lb[1], node$lb[1]), y = c(node$ub[2], node$lb[2]))
    points(node$centroid[[1]], node$centroid[[2]], col = col, pch = 20, cex = 0.5)
    # plotrix::draw.circle(x = node$centroid[[1]], y = node$centroid[[2]], radius = node$radius)
    text(x = node$centroid[[1]], y = node$centroid[[2]], labels = node$id, col = col, pos = 3, cex = 0.5)
  }
}

plot(x, xlim = range(x[, 1]) + c(-1, 1), ylim = range(x[, 2]) + c(-1, 1), asp = 1)
text(x, labels = 1:nrow(x), pos = 3)
plotTreeProp(info$tree_properties)

microbenchmark::microbenchmark({ info <- clustertree:::testMST(x2, bucketSize = 30L) }, times = 1L)
microbenchmark::microbenchmark({ wut <- dbscan:::prims(dist(x2), nrow(x2)) }, times = 1L)

  info <- clustertree:::testMST(x, bucketSize = 5L)
  (sum(dbscan:::prims(dist(x), nrow(x))[, 3]) - sum(sqrt(info$mst[, 3]))) < sqrt(.Machine$double.eps)

## 2.27183
## Testing tree construction time
  clustertree:::testTrees(x2, k = 25, bucketSize = 10) ## 11,19 should be pruned

## Profiling
  microbenchmark::microbenchmark({ dbscan::kNN(x2, k = 25, bucketSize = 10) }, times = 1L)
  microbenchmark::microbenchmark({ clustertree:::testDTS(x2, k = 25, bucketSize = 10) }, times = 1L)

  invisible(gprofiler::profile({ replicate(5, invisible(clustertree:::testDTS(x2, k = 25, bucketSize = 10))) }, filename = "cluster_perf.out"))
  system(command = "pprof.pl --web /Library/Frameworks/R.framework/Resources/bin/exec/R /Users/mpiekenbrock/clustertree/clustertree_prof.out")

  plotLeafBB <- function(nodes, leaf_only = FALSE, col = "red"){
    for (node in nodes){
      lines(x = c(node$lb[1], node$ub[1]), y = c(node$lb[2], node$lb[2]))
      lines(x = c(node$ub[1], node$ub[1]), y = c(node$lb[2], node$ub[2]))
      lines(x = c(node$ub[1], node$lb[1]), y = c(node$ub[2], node$ub[2]))
      lines(x = c(node$lb[1], node$lb[1]), y = c(node$ub[2], node$lb[2]))
      points(node$centroid[[1]], node$centroid[[2]], col = col, pch = 20, cex = 0.5)
      plotrix::draw.circle(x = node$centroid[[1]], y = node$centroid[[2]], radius = node$radius)
    }
  }
plot(x, xlim = range(x[, 1]) + c(-1, 1), ylim = range(x[, 2]) + c(-1, 1), asp = 1)
  text(x, labels = 0:(nrow(x) - 1), pos = 3)

  node_ids <- sapply(info$info, function(node) node$id)
  leaf_idx <- sapply(info$info, function(node) !is.null(node$idx))
  plotLeafBB(info$info[leaf_idx])

  split_idx <- sapply(info$info, function(node) is.null(node$idx))
  split_centroids <- t(sapply(info$info[split_idx], function(node) node$centroid))
  split_id_int <- c(0)
  split_int <- which(node_ids[split_idx] %in% split_id_int)
  plotLeafBB(info$info[which(node_ids %in% split_id_int)], col = "blue")
  text(split_centroids[split_int,1], split_centroids[split_int,2], labels = node_ids[split_idx][split_int], pos = 3, col = "blue", cex = 0.65)

  leaf_centroids <- t(sapply(info$info[leaf_idx], function(node) node$centroid))
  text(leaf_centroids, labels = node_idx[leaf_idx], pos = 3, col = "red", cex=.5)

  which(sapply(sapply(info$info, function(node) node$id), function(id) if(is.null(id)) FALSE else id == 1L))
  score(info$info[[8]], info$info[[2]])


  score <- function(qn, rn){
    dist_qr <- dist(rbind(qn$centroid, rn$centroid))[[1]]^2
    dist_qr - qn$max_radius - rn$max_radius
  }

xb <- cbind(rnorm(1000), rnorm(1000)) #as.matrix(iris[, 1:4]) #
  cl1 <- function(){ invisible(clustertree:::testDTS(xb, k = 15, bucketSize = 15L)) }
cl2 <- function(){ invisible(dbscan::kNN(xb, k = 15, bucketSize = 15L)) }
microbenchmark::microbenchmark(cl1(), times = 3L)
  microbenchmark::microbenchmark(cl2(), times = 3L)

  all((cl2()$idx - 1) == cl1()$idx)

  */
