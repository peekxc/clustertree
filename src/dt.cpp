#include <Rcpp.h>
using namespace Rcpp;

// Local header includes
#include "dt.h"

// DualTree enables the dual-traversal of two kd-tree's at the same time, and
// for computing the associated bounds that come with them
DualTree::DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree) : d(ref_tree->theDim()) {
  if (ref_tree->theDim() != query_tree->theDim()){ stop("Dimensionality of the query set does not match the reference set."); }
  rtree = ref_tree, qtree = query_tree;  // Store pointers to both trees
}

// TODO: Export better to R, figure out how to extract R output buffer instead of cout for windows users
void DualTree::PrintTree(bool ref_tree){
  if (ref_tree){ rtree->Print((ANNbool) true, std::cout); }
  else { qtree->Print((ANNbool) true, std::cout); }
}


IntegerVector DualTree::getIDXArray(){
  IntegerVector ids = Rcpp::no_init(rtree->n_pts);
  for (int i = 0; i < rtree->n_pts; ++i){
    ids[i] = rtree->pidx[i];
  }
  return ids;
}

// The minimum distance between any two points descendent of N_i and N_j
ANNdist DualTree::min_dist(ANNkd_node* N_i, ANNkd_node* N_j){
  if (N_i == N_j) return (ANNdist) 0.0; // Trivial case--nodes are the same
  //Rcout << "Computing the minimum distance" << std::endl;
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->desc_ids(*i_desc_ids);
  N_j->desc_ids(*j_desc_ids);

  ANNdist min_dist = ANN_DIST_INF, dist;
  for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
    for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
      ANNpoint p_i = this->rtree->pts[*i];
      ANNpoint p_j = this->qtree->pts[*j];
      dist = annDist(d, p_i, p_j);
      if (dist < min_dist){
        min_dist = dist;
      }
    }
  }
  return min_dist;
}

// The maximum distance between any two points descendent of N_i and N_j
ANNdist DualTree::max_dist(ANNkd_node* N_i, ANNkd_node* N_j){
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->desc_ids(*i_desc_ids);
  N_j->desc_ids(*j_desc_ids);

  ANNdist max_dist = 0, dist;
  for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
    for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
      ANNpoint p_i = this->rtree->pts[*i];
      ANNpoint p_j = this->qtree->pts[*j];
      if ((dist = annDist(d, p_i, p_j)) > max_dist){
        max_dist = dist;
      }
    }
  }
  return max_dist;
}

// Get the maximum distance between the centroid of the convex subset (the box of the current node)
// and the points within the children of the current node
ANNdist DualTree::max_child_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){

  // If the maximum child distance has already been computed, return it
  if (ni_bnd.rho != -1.0) return ni_bnd.rho;

  // Else look at the points its holding
  std::vector<int> point_ids = std::vector<int>();
  N_i->node_ids(point_ids);

  // Centroid computes (and stores) the bounding box, so access should be fine past here
  ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;

  // If it's a splitting node it contains no points. There are two options:
  // 1) Search all of the descendents, computing the true maximum child distance
  // 2) Approximate a looser upper bound using the distance between the precomputed centroid
  // and the bounding rectangle orthogonal to the cutting dimension
  if (point_ids.size() == 0){
    ni_bnd.rho = annDist(d, ni_centroid, ni_bnd.bnd_box->hi); // <-- Looser upper bound
    // Rcout << "== Computed child distance: " << ni_bnd.rho << std::endl;
    // ni_bnd.rho = max_desc_dist(N_i, ni_bnd, ref_tree); // <-- Tighter but slower to compute bound
    return ni_bnd.rho;
  } else {
    // If it's a leaf node and it contains no points, find the one farthest from the centroid
    ANNdist max_dist = 0;
    for (std::vector<int>::iterator pt_idx = point_ids.begin(); pt_idx != point_ids.end(); ++pt_idx){
      ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
      ANNdist dist = annDist(d, cpt, ni_centroid);
      // Rcout << "Testing child distance: " << dist << std::endl;
      if (dist > max_dist){ max_dist = dist; }
    }
    ni_bnd.rho = max_dist; // store the bound
    return max_dist;
  }
}

// Get the maximum distance between the centroid of the convex subset (the box of the current node)
// and all of the descendent points within a given node
ANNdist DualTree::max_desc_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){

  // If the maximum desc. distance has already been computed, return it
  if (ni_bnd.lambda != -1.0) return ni_bnd.lambda;

  // Get the centroid
  ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;

  // Get all descendent point indices
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNdist max_dist = 0;
  bool check_point_knn = false;
  for (std::vector<int>::iterator pt_idx = desc_ids.begin(); pt_idx != desc_ids.end(); ++pt_idx){
    ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
    ANNdist dist = annDist(d, cpt, ni_centroid);
    if (dist > max_dist){ max_dist = dist; }
    // if (knn->at(*pt_idx)->max_key() != ANN_DIST_INF){
    //   check_point_knn = true;
    //   Rcout << "*** " << N_i << " has finite knn points filled!" << std::endl;
    // }
  }
  ni_bnd.lambda = max_dist; // Store the maximum distance by reference for later
  return max_dist;
}


#include "kd_util.h"
#include "ANN/ANNx.h"
#include "ANN/ANNperf.h"

ANNorthRect& DualTree::convex_subset(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
  if (ni_bnd.bnd_box != NULL) return *ni_bnd.bnd_box;
  // Get child ids
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNidxArray pidx = (ANNidxArray) &desc_ids[0];

  // Compute the bounding rectangle
  ANNorthRect* bnd_box = new ANNorthRect(d); // bounds are dynamically allocated
  annEnclRect((ANNpointArray) ref_tree ? rtree->pts : qtree->pts,
              (ANNidxArray) pidx,
              (int) desc_ids.size(),
              (int) d,
              *bnd_box);
  ni_bnd.bnd_box = bnd_box; // Store bounding rectangle
  return *bnd_box;
}

ANNpoint DualTree::centroid(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
  if (ni_bnd.centroid != NULL) return(ni_bnd.centroid);
  ANNorthRect& bnds = ni_bnd.bnd_box == NULL ? convex_subset(N_i, ni_bnd, ref_tree) : *ni_bnd.bnd_box;
  ANNpoint n_centroid = new ANNcoord[d];
  for (int i = 0; i < d; ++i){
    // Rcout << "Bounds (lo, hi): " << bnds.lo[i] << ", " << bnds.hi[i] << std::endl;
    n_centroid[i] = ANNcoord((bnds.lo[i] + bnds.hi[i]) / 2.0);
    // Rcout << "Centroid: " << centroid[i] << std::endl;
  }
  ni_bnd.centroid = n_centroid;
  return(n_centroid);
}

/*** R

  plot(test_set, xlim = range(test_set[, 1]) + c(-50, 50), ylim = range(test_set[, 2]) + c(-50, 50),
       asp = 1)
  text(test_set, labels = 1:10, pos = 3)
  clustertree:::DT_knn(test_set, k = 4L, bkt_size = 2L)$id[, -1] == unname(dbscan::kNN(test_set, k = 4L)$id)


  clustertree:::DT_knn(test_set, k = 3L)$dist[, -1] == unname(dbscan::kNN(test_set, k = 3L)$dist)



  size <- 1500
  xyz <- as.matrix(data.frame(x = rnorm(size), y = rnorm(size), z = rnorm(size)))
  microbenchmark::microbenchmark(invisible(FNN::get.knn(xyz, k = 30L, algorithm = "kd_tree")), times = 15L)
  microbenchmark::microbenchmark(invisible(clustertree:::DT_knn(xyz, k = 30L, bkt_size = 15L, prune = TRUE)), times = 15L)
  microbenchmark::microbenchmark(invisible(clustertree:::DT_knn(xyz, k = 30L, bkt_size = 15L, prune = FALSE)), times = 15L)
  microbenchmark::microbenchmark(invisible(dbscan::kNN(xyz, k = 30L, sort = F, bucketSize = 15L)), times = 15L)

  */


