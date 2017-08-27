/*  Unexported, experimental functions.
 */

#include <Rcpp.h>
using namespace Rcpp;

#include <ANN/ANN_util.h> // matrixToANNpointArray
#include <clustertree/greedy_prims.h>

List dist_prims(const NumericMatrix& x) {

  // Copy data over to ANN point array
  int n = x.nrow(), d = x.ncol();
  ANNpointArray dataPts = matrixToANNpointArray(x);

  // Create kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, n, d, 10, ANN_KD_SUGGEST);

  // Distances
  NumericVector dist = no_init(n - 1);
  IntegerMatrix merge = no_init_matrix(n - 1, 2);

  // Data structures for prims
  // The 'fringe' maps the shortest known distance of the 'from' node by index, to the 'to' node
  std::vector<int> visited = std::vector<int>(n, -1); // -1 to indicate node is not in MST
  std::vector<simple_edge> fringe = std::vector<simple_edge>(n, simple_edge(-1, ANN_DIST_INF));

  // Current mind distance and min id
  double min;
  int min_id = 0;

  // Start with the current node, which should be unvisited; Find it's nearest unvisited NN. Update the fringe for it.
  simple_edge se;// = findNearest(c_i, kdTree, visited);
  visited[min_id] = 1;
  for (int n_edges = 0; n_edges < n - 1; n_edges++) {
    if (n_edges % 1000 == 0) Rcpp::checkUserInterrupt();
    min = ANN_DIST_INF; // Reset every iteration to infinity so new edge is always chosen

    // Need to now update all of the nodes that have been visited before
    for (int v_i = 0; v_i < n; ++v_i) {
      // Update previously visited node NN's, excluding the current candidate point!
      if (visited[v_i] >= 0) {
        se = findNearest(v_i, kdTree, visited);
        fringe[v_i] = se; // this is the lowest weight for v_i
      }

      // Prim's step: determine whether to continue along this node or a node on the fringe
      int to_node = fringe[v_i].to;
      if (fringe[v_i].weight < min && visited[to_node] < 0){ // every visited node should already have the updated non-visited nearest
        min = fringe[v_i].weight, min_id = v_i; // min_id is the node containing the shortest edge to a non-visited node
      }
    }

    // Record the from id (min_id) and it's target minimum (fringe[min_id].to)
    merge(n_edges, _) = IntegerVector::create(min_id, fringe[min_id].to);
    dist[n_edges] = min;
    visited[min_id] = visited[fringe[min_id].to] = 1;
  }

  // Cleanup
  delete kdTree;
  // delete [] dists;
  // delete [] nnIdx;
  annDeallocPts(dataPts);
  annClose();

  return(List::create(_["merge"] = merge, _["dist"] = dist));
}


List knn_test(const NumericMatrix& x, const int k) {

  // Copy data over to ANN point array
  int n = x.nrow(), d = x.ncol();
  ANNpointArray dataPts = matrixToANNpointArray(x);

  // Create kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, n, d, 10, ANN_KD_SUGGEST);

  NumericMatrix dist(n, k);
  IntegerMatrix id(n, k);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  for (int i=0; i < n; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = dataPts[i];
    kdTree->annkSearch(queryPt, k+1, nnIdx, dists, 0.0);

    NumericVector ndists = NumericVector(dists, dists+k+1);
    dist(i, _) = sqrt(ndists);
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k+1);
    id(i, _) = ids + 1;
  }

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(dataPts);
  annClose();

  // prepare results
  List ret;
  ret["dist"] = d;
  ret["id"] = id;
  ret["k"] = k;
  return ret;
}

List knn_test2(const NumericMatrix& x, const int k) {

  // Copy data over to ANN point array
  int n = x.nrow(), d = x.ncol();
  ANNpointArray dataPts = matrixToANNpointArray(x);

  // Create kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, n, d, 10, ANN_KD_SUGGEST);

  NumericMatrix dist(n, k);
  IntegerMatrix id(n, k);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  for (int i=0; i<n; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();
    knn_stack(i, k,	nnIdx, dists,	kdTree);
    NumericVector ndists = NumericVector(dists, dists+k); // don't add the plus 1
    dist(i, _) = sqrt(ndists);
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k);// don't add the plus 1
    id(i, _) = ids + 1;
  }

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(dataPts);
  annClose();

  // prepare results
  List ret;
  ret["dist"] = d;
  ret["id"] = id;
  ret["k"] = k;
  return ret;
}


